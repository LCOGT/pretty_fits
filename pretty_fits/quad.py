import numpy as np
import os
from astropy.io import fits
import operator
import itertools
import logging
from scipy.spatial.distance import cdist
import math

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def calibrate(refimg, imgs):

    logger.debug("Analysing reference image")
    ref = ImgCat(refimg)
    ref.makestarlist()
    ref.makemorequads()

    identifications = []

    logger.debug("Analysing images from stack")
    for ukn in imgs:

        ukn = ImgCat(ukn)
        ukn.makestarlist()

        idn = Identification(ref, ukn)
        idn.findtrans()
        idn.calcfluxratio()
        identifications.append(idn)

    return identifications

class Identification:
    """
    Represents the identification of a transform between two ImgCat objects.
    Regroups all the star catalogs, the transform, the quads, the candidate, etc.

    All instance attributes are listed below.

    :ivar ref: ImgCat object of the reference image
    :ivar ukn: ImgCat object of the unknown image
    :ivar ok: boolean, True if the idendification was successful.
    :ivar trans: The SimpleTransform object that represents the geometrical transform from ukn to ref.
    :ivar uknmatchstars: A list of Star objects of the catalog of the unknown image...
    :ivar refmatchstars: ... that correspond to these Star objects of the reference image.
    :ivar medfluxratio: Median flux ratio between the images: "ukn * medfluxratio = ref"
        A high value corresponds to a shallow unknown image.
        It gets computed by the method calcfluxratio, using the matched stars.
    :ivar stdfluxratio: Standard error on the flux ratios of the matched stars.



    """

    def __init__(self, ref, ukn):
        """

        :param ref: The reference image
        :type ref: ImgCat object
        :param ukn: The unknown image, whose transform will be adjusted to match the ref
        :type ukn: ImgCat object

        """
        self.ref = ref
        self.ukn = ukn

        self.ok = False

        self.trans = None
        self.uknmatchstars = []
        self.refmatchstars = []
        self.cand = None

        self.medfluxratio = None # ukn * medfluxratio --> ref (high value means shallow image)
        self.stdfluxratio = None


    def findtrans(self, r = 5.0, verbose=True):
        """
        Find the best trans given the quads, and tests if the match is sufficient
        """

        # Some robustness checks
        if len(self.ref.starlist) < 4:
            logger.debug("Not enough stars in the reference catalog.")
            return
        if len(self.ukn.starlist) < 4:
            logger.debug("Not enough stars in the unknown catalog.")
            return


        # First question : how many stars should match ?
        if len(self.ukn.starlist) < 5: # Then we should simply try to get the smallest distance...
            minnident = 4
        else:
            minnident = max(4, min(8, len(self.ukn.starlist)/5.0)) # Perfectly arbitrary, let's see how it works

        # Hmm, arbitrary for now :
        minquaddist = 0.005

        # Let's start :
        if self.ref.quadlevel == 0:
            self.ref.makemorequads(verbose=verbose)
        if self.ukn.quadlevel == 0:
            self.ukn.makemorequads(verbose=verbose)

        while self.ok == False:

            # Find the best candidates
            cands = proposecands(self.ukn.quadlist, self.ref.quadlist, n=4, verbose=verbose)

            if len(cands) != 0 and cands[0]["dist"] < minquaddist:
                # If no quads are available, we directly try to make more ones.
                for cand in cands:
                    # Check how many stars are identified...
                    nident = identify(self.ukn.starlist, self.ref.starlist, trans=cand["trans"], r=r, verbose=verbose, getstars=False)
                    if nident >= minnident:
                        self.trans = cand["trans"]
                        self.cand = cand
                        self.ok = True
                        break # get out of the for

            if self.ok == False:
                # We add more quads...
                addedmorerefquads = self.ref.makemorequads(verbose=verbose)
                addedmoreuknquads = self.ukn.makemorequads(verbose=verbose)

                if addedmorerefquads == False and addedmoreuknquads == False:
                    break # get out of the while, we failed.


        if self.ok: # we refine the transform
            # get matching stars :
            (self.uknmatchstars, self.refmatchstars) = identify(self.ukn.starlist, self.ref.starlist, trans=self.trans, r=r, verbose=False, getstars=True)
            # refit the transform on them :
            logger.debug("Refitting transform (before/after) :")
            logger.debug(self.trans)
            newtrans = fitstars(self.uknmatchstars, self.refmatchstars)
            if newtrans != None:
                self.trans = newtrans
                logger.debug(self.trans)
            # Generating final matched star lists :
            (self.uknmatchstars, self.refmatchstars) = identify(self.ukn.starlist, self.ref.starlist, trans=self.trans, r=r, verbose=verbose, getstars=True)

            logger.debug("I'm done !")

        else:
            logger.debug("Failed to find transform !")

    def calcfluxratio(self, verbose=True):
        """
        Computes a very simple median flux ratio between the images.
        The purpose is simply to get a crude guess, for images with e.g. different exposure times.
        Given that we have these corresponding star lists in hand, this is trivial to do once findtrans was run.
        """
        assert len(self.uknmatchstars) == len(self.refmatchstars)
        if len(self.refmatchstars) == 0:
            logger.debug("No matching stars to compute flux ratio !")
            return

        reffluxes = listtoarray(self.refmatchstars, full=True)[:,2]
        uknfluxes = listtoarray(self.uknmatchstars, full=True)[:,2]
        fluxratios = reffluxes / uknfluxes

        self.medfluxratio = float(np.median(fluxratios))
        self.stdfluxratio = float(np.std(fluxratios))

        logger.debug("Computed flux ratio from %i matches : median %.2f, std %.2f" % (len(reffluxes), self.medfluxratio, self.stdfluxratio))


class ImgCat:
    """
    Represent an individual image and its associated catalog, starlist, quads etc.
    """

    def __init__(self, filepath, hdu=0, cat=None):
        """

        :param filepath: Path to the FITS file, or alternatively just a string to identify the image.
        :type filepath: string

        :param cat: Catalog generated by SExtractor (if available -- if not, we'll make our own)
        :type cat: asciidata catalog

        :param hdu: The hdu containing the science data from which I should build the catalog. 0 is primary. If multihdu, 1 is usually science.

        """
        self.filepath = filepath

        (imgdir, filename) = os.path.split(filepath)
        (common, ext) = os.path.splitext(filename)
        self.name = common

        self.hdu = hdu
        self.cat = cat
        self.starlist = []
        self.mindist = 0.0
        self.xlim = (0.0, 0.0) # Will be set using the catalog -- no need for the FITS image.
        self.ylim = (0.0, 0.0)

        self.quadlist = []
        self.quadlevel = 0 # encodes what kind of quads have already been computed

    def makestarlist(self, skipsaturated=False, n=200):
        if skipsaturated:
            maxflag = 3
        else:
            maxflag = 7
        hdu = fits.open(self.filepath)
        cats = hdu[2].data
        self.starlist = sortstarlistbyflux(cats)[:n]
        (xmin, xmax, ymin, ymax) = area(cats, border=0.01)
        self.xlim = (xmin, xmax)
        self.ylim = (ymin, ymax)

        # Given this starlists, what is a good minimal distance for stars in quads ?
        self.mindist = min(min(xmax - xmin, ymax - ymin) / 10.0, 30.0)

    def makemorequads(self, verbose=True):
        """
        We add more quads, following the quadlevel.
        """
        #if not add:
        #    self.quadlist = []
        if verbose:
            logger.debug("Making more quads, from quadlevel %i ..." % self.quadlevel)
        if self.quadlevel == 0:
            self.quadlist.extend(makequads1(self.starlist, n=7, d=self.mindist, verbose=verbose))
        elif self.quadlevel == 1:
            self.quadlist.extend(makequads2(self.starlist, f=3, n=5, d=self.mindist, verbose=verbose))
        elif self.quadlevel == 2:
            self.quadlist.extend(makequads2(self.starlist, f=6, n=5, d=self.mindist, verbose=verbose))
        elif self.quadlevel == 3:
            self.quadlist.extend(makequads2(self.starlist, f=12, n=5, d=self.mindist, verbose=verbose))
        elif self.quadlevel == 4:
            self.quadlist.extend(makequads2(self.starlist, f=10, n=6, s=3, d=self.mindist, verbose=verbose))

        else:
            return False

        self.quadlist = removeduplicates(self.quadlist, verbose=verbose)
        self.quadlevel += 1
        return True

class Quad:
    """
    A geometric "hash", or asterism, as used in Astrometry.net :
    http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0910.2233
    It is made out of 4 stars, and it is shift / scale / rotation invariant
    """

    def __init__(self, fourstars):
        """
        fourstars is a list of four stars

        We make the following attributes :
        self.hash
        self.stars (in the order A, B, C, D)

        """
        assert len(fourstars) == 4

        tests = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
        other = [(2,3), (1,3), (1,2), (0,3), (0,2), (0,1)]
        dists = np.array([np.linalg.norm(np.array(fourstars[0]['x'], fourstars[0]['y']) - np.array(fourstars[1]['x'], fourstars[1]['y'])) for (i,j) in tests])
        assert np.min(dists) > 1.0

        maxindex = np.argmax(dists)
        (Ai, Bi) = tests[maxindex] # Indexes of stars A and B
        (Ci, Di) = other[maxindex] # Indexes of stars C and D
        A = fourstars[Ai]
        B = fourstars[Bi]
        C = fourstars[Ci]
        D = fourstars[Di]

        # We look for matrix transform [[a -b], [b a]] + [c d] that brings A and B to 00 11 :
        x = B['x'] - A['x']
        y = B['y'] - A['y']
        b = (x-y)/(x*x + y*y)
        a = (1.0/x) * (1.0 + b*y)
        c = b*A['y'] - a*A['x']
        d = - (b*A['x'] + a*A['y'])

        t = SimpleTransform((a, b, c, d))

        # Test
        #logger.debug(t.apply((A['x'], A['y'])))
        #logger.debug(t.apply((B.x, B['y'])))

        (xC, yC) = t.apply(x = C['x'], y = C['y'])
        (xD, yD) = t.apply(x = D['x'], y = D['y'])

        # Normal case
        self.hash = (xC, yC, xD, yD)

        # Break symmetries :
        testa = xC > xD
        testb = xC + xD > 1

        if testa and not testb: # we switch C and D
            #logger.debug("a")
            self.hash = (xD, yD, xC, yC)
            (C, D) = (D, C)

        if testb and not testa: # We switch A and B
            #logger.debug("b")
            self.hash = (1.0-xD, 1.0-yD, 1.0-xC, 1.0-yC)
            (A, B) = (B, A)
            (C, D) = (D, C)

        if testa and testb:
            #logger.debug("a + b")
            self.hash = (1.0-xC, 1.0-yC, 1.0-xD, 1.0-yD)
            (A, B) = (B, A)

        # Checks :
        assert self.hash[0] <= self.hash[2]
        assert self.hash[0] + self.hash[2] <= 1

        self.stars = [A, B, C, D] # Order might be different from the fourstars !


    def __str__(self):
        return "Hash : %6.3f %6.3f %6.3f %6.3f / IDs : (%s, %s, %s, %s)" % (
            self.hash[0], self.hash[1], self.hash[2], self.hash[3],
            self.stars[0].name, self.stars[1].name, self.stars[2].name, self.stars[3].name)

class SimpleTransform:
    """
    Represents an affine transformation consisting of rotation, isotropic scaling, and shift.
    [x', y'] = [[a -b], [b a]] * [x, y] + [c d]
    """

    def __init__(self, v = (1, 0, 0, 0)):
        """
        v = (a, b, c, d)
        """
        self.v = np.asarray(v)

    def getscaling(self):
        return math.sqrt(self.v[0]*self.v[0] + self.v[1]*self.v[1])

    def getrotation(self):
        """
        The CCW rotation angle, in degrees
        """
        return math.atan2(self.v[1], self.v[0]) * (180.0/math.pi)# % 360.0

    def __str__(self):
        return "Rotation %+11.6f [deg], scale %8.6f" % (self.getrotation(), self.getscaling())


    def inverse(self):
        """
        Returns the inverse transform !
        """

        # To represent affine transformations with matrices, we can use homogeneous coordinates.
        homo = np.array([
        [self.v[0], -self.v[1], self.v[2]],
        [self.v[1],  self.v[0], self.v[3]],
        [0.0, 0.0, 1.0]
        ])

        inv = np.linalg.inv(homo)
        #logger.debug(inv)

        return SimpleTransform((inv[0,0], inv[1,0], inv[0,2], inv[1,2]))



    def matrixform(self):
        """
        Special output for scipy.ndimage.interpolation.affine_transform
        Returns (matrix, offset)
        """

        return (np.array([[self.v[0], -self.v[1]], [self.v[1], self.v[0]]]), self.v[2:4])


    def apply(self, x, y):
        """
        Applies the transform to a point (x, y)
        """
        xn = self.v[0]*x -self.v[1]*y + self.v[2]
        yn = self.v[1]*x +self.v[0]*y + self.v[3]
        return (xn, yn)

    def applystar(self, transstar):
        (transstar['x'], transstar['y']) = self.apply(x=transstar['x'], y=transstar['y'])
        return transstar

    def applystarlist(self, starlist):
        return [self.applystar(star) for star in starlist]

def sortstarlistbyflux(starlist):
    """
    We sort starlist according to flux : highest flux first !
    """
    sortedstarlist = sorted(starlist, key=operator.itemgetter('flux'))
    sortedstarlist.reverse()
    return sortedstarlist

def area(stars, border=0.01):
    """
    Returns the area covered by the stars.
    Border is relative to max-min
    """
    if len(stars) == 0:
        return np.array([0, 1, 0, 1])

    if len(stars) == 1:
        star = stars[0]
        return np.array([star[0] - 0.5, star[0] + 0.5, star[1] - 0.5, star[1] + 0.5])

    (xmin, xmax) = (np.min(stars[0]), np.max(stars[0]))
    (ymin, ymax) = (np.min(stars[1]), np.max(stars[1]))
    xw = xmax - xmin
    yw = ymax - ymin
    xmin = xmin - border*xw
    xmax = xmax + border*xw
    ymin = ymin - border*yw
    ymax = ymax + border*yw
    return np.array([xmin, xmax, ymin, ymax])

def makequads1(starlist, n=7, s=0, d=50.0, verbose=True):
    """
    First trivial quad maker.
    Makes combis of the n brightest stars.

    :param n: number of stars to consider (brightest ones).
    :type n: int
    :param s: how many of the brightest stars should I skip ?
        This feature is useful to avoid building quads with nearly saturated stars that are not
        available in other exposures.
    :type s: int
    :param d: minimal distance between stars
    :type d: float

    """
    quadlist = []
    sortedstars = sortstarlistbyflux(starlist)

    for fourstars in itertools.combinations(sortedstars[s:s+n], 4):
        if mindist(fourstars) > d:
                quadlist.append(Quad(fourstars))

    logger.debug("Made %4i quads from %4i stars (combi n=%i s=%i d=%.1f)" % (len(quadlist), len(starlist), n, s, d))

    return quadlist

def makequads2(starlist, f=5.0, n=6, s=0, d=50.0, verbose=True):
    """
    Similar, but fxf in subareas roughly f times smaller than the full frame.
    s allows to skip the brightest stars in each region

    :param f: smallness of the subareas
    :type f: float
    :param n: number of stars to consider in each subarea
    :type n: int
    :param d: minimal distance between stars
    :type d: float
    :param s: number of brightest stars to skip in each subarea
    :type s: int

    """
    quadlist = []
    sortedstars = sortstarlistbyflux(starlist)
    (xmin, xmax, ymin, ymax) = area(sortedstars)

    r = 2.0 * max(xmax - xmin, ymax - ymin) / f

    for xc in np.linspace(xmin, xmax, f+2)[1:-1]:
        for yc in np.linspace(ymin, ymax, f+2)[1:-1]:
            cstar = (xc, yc)
            das = distanceandsort(cstar, sortedstars)
            #closest = [s["star"] for s in das[0:4]]
            brightestwithinr = sortstarlistbyflux([element["star"] for element in das if element["dist"] <= r])[s:s+n]
            for fourstars in itertools.combinations(brightestwithinr, 4):
                if mindist(fourstars) > d:
                    quadlist.append(Quad(fourstars))

    logger.debug("Made %4i quads from %4i stars (combi sub f=%.1f n=%i s=%i d=%.1f)" % (len(quadlist), len(starlist), f, n, s, d))

    return quadlist

def mindist(cats):
    """
    Function that tests if 4 stars are suitable to make a good quad...
    """
    tests = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
    dists = np.array([np.linalg.norm(np.array(cats[0]['x'], cats[0]['y']) - np.array(cats[1]['x'], cats[1]['y'])) for (i,j) in tests])
    return np.min(dists)

def removeduplicates(quadlist, verbose=True):
    """
    Returns a quadlist without quads with identical hashes...
    """
    # To avoid crash in lexsort if quadlist is too small :
    if len(quadlist) < 2:
        return quadlist
    hasharray = np.array([q.hash for q in quadlist])

    order = np.lexsort(hasharray.T)
    hasharray = hasharray[order]
    #diff = np.diff(hasharray, axis=0)
    diff = np.fabs(np.diff(hasharray, axis=0))
    #diff = np.sum(diff, axis=1)
    ui = np.ones(len(hasharray), 'bool')
    ui[1:] = (diff >= 0.000001).any(axis=1)
    #logger.debug(hasharray[ui==False])
    logger.debug("Removing %i/%i duplicates" % (len(quadlist) - np.sum(ui), len(quadlist)))

    return [quad for (quad, u) in zip(quadlist, ui) if u == True]

def proposecands(uknquadlist, refquadlist, n=5, verbose=True):
    """
    Function that identifies similar quads between the unknown image and a reference.
    Returns a dict of (uknquad, refquad, dist, trans)
    """
    # Nothing to do if the quadlists are empty ...
    if len(uknquadlist) == 0 or len(refquadlist) == 0:
        logger.debug("No quads to propose ...")
        return []

    logger.debug("Finding %i best candidates among %i x %i (ukn x ref)" % (n, len(uknquadlist), len(refquadlist)))
    uknhashs = np.array([q.hash for q in uknquadlist])
    refhashs = np.array([q.hash for q in refquadlist])

    # Brute force...
    dists = cdist(refhashs, uknhashs)
    uknmindistindexes = np.argmin(dists, axis=0) # For each ukn, the index of the closest ref
    uknmindist = np.min(dists, axis=0) # The corresponding distances
    uknbestindexes = np.argsort(uknmindist)

    candlist = []
    nmax = len(uknbestindexes)
    logger.debug("We have a maximum of %i quad pairs" % (nmax))
    logger.debug("{} {}".format(uknbestindexes, dists))
    for i in range(min(n, nmax)):

        cand = {"uknquad": uknquadlist[uknbestindexes[i]], "refquad":refquadlist[uknmindistindexes[uknbestindexes[i]]],
            "dist":uknmindist[uknbestindexes[i]]}

        cand["trans"] = quadtrans(cand["uknquad"], cand["refquad"])

        candlist.append(cand)
        logger.debug("Cand %2i (dist. %12.8f) : %s" % (i+1, cand["dist"], str(cand["trans"])))

    return candlist

def fitstars(uknstars, refstars, verbose=True):
    """
    I return the transform that puts the unknown stars (uknstars) onto the refstars.
    If you supply only two stars, this is using linalg.solve() -- perfect solution.
    If you supply more stars, we use linear least squares, i.e. minimize the 2D error.

    Formalism inspired by :
    http://math.stackexchange.com/questions/77462/
    """

    assert len(uknstars) == len(refstars)
    if len(uknstars) < 2:
        logger.debug("Sorry I cannot fit a transform on less than 2 stars.")
        return None

    # ukn * x = ref
    # x is the transform (a, b, c, d)

    ref = np.hstack(listtoarray(refstars)) # a 1D vector of lenth 2n

    uknlist = []
    for star in uknstars:
        uknlist.append([star['x'], -star['y'], 1, 0])
        uknlist.append([star['y'], star['x'], 0, 1])
    ukn = np.vstack(np.array(uknlist)) # a matrix

    if len(uknstars) == 2:
        trans = np.linalg.solve(ukn, ref)
    else:
        trans = np.linalg.lstsq(ukn, ref)[0]

    return SimpleTransform(np.asarray(trans))

def quadtrans(uknquad, refquad):
    """
    Quickly return a transform estimated from the stars A and B of two quads.
    """
    return fitstars(uknquad.stars[:2], refquad.stars[:2])

def listtoarray(starlist, full=False):
    """
    Transforms the starlist into a 2D numpy array for fast manipulations.
    First index is star, second index is x or y

    :param full: If True, I also include flux, fwhm, ELLIPTICITY
    :type full: boolean

    """
    if full:
        return np.array([[star['x'], star['y'], star['flux'], star['fwhm'], star['ELLIPTICITY']] for star in starlist])
    else:
        return np.array([[star['x'], star['y']] for star in starlist])

def identify(uknstars, refstars, trans=None, r=5.0, verbose=True, getstars=False):
    """
    Allows to:
     * get the number or matches, i.e. evaluate the quality of the trans
     * get corresponding stars from both lists (without the transform applied)

    :param getstars: If True, I return two lists of corresponding stars, instead of just the number of matching stars
    :type getstars: boolean

    Inspired by the "formpairs" of alipy 1.0 ...

    """

    if trans != None:
        ukn = listtoarray(trans.applystarlist(uknstars))
    else:
        ukn = listtoarray(uknstars)
    ref = listtoarray(refstars)

    dists = cdist(ukn, ref) # Big table of distances between ukn and ref
    mindists = np.min(dists, axis=1) # For each ukn, the minimal distance
    minok = mindists <= r # booleans for each ukn
    minokindexes = np.argwhere(minok).flatten() # indexes of uknstars with matches

    if verbose:
        logger.debug("%i/%i stars with distance < r = %.1f (mean %.1f, median %.1f, std %.1f)" % (np.sum(minok), len(uknstars), r,
            np.mean(mindists[minok]), np.median(mindists[minok]), np.std(mindists[minok])))

    matchuknstars = []
    matchrefstars = []

    for i in minokindexes: # we look for the second nearest ...
        sortedrefs = np.argsort(dists[i,:])
        firstdist = dists[i,sortedrefs[0]]
        seconddist = dists[i,sortedrefs[1]]
        if seconddist > 2.0*firstdist: # Then the situation is clear, we keep it.
            matchuknstars.append(uknstars[i])
            matchrefstars.append(refstars[sortedrefs[0]])
        else:
            pass # Then there is a companion, we skip it.

    if verbose:
        logger.debug("Filtered for companions, keeping %i/%i matches" % (len(matchuknstars), np.sum(minok)))

    if getstars==True:
        return (matchuknstars, matchrefstars)
    else:
        return len(matchuknstars)

def distanceandsort(cstar, otherstarlist):
    """
    Returns a list of dicts(star, dist, origpos), sorted by distance to self.
    The 0th star is the closest.

    otherstarlist is not modified.
    """
    import operator # for the sorting

    returnlist=[]
    for i, star in enumerate(otherstarlist):
        dist = np.sqrt(np.sum((np.array(cstar) - np.array([star[0],star[1]]))**2))#self.distance(star)
        returnlist.append({'star':star, 'dist':dist, 'origpos':i})
    returnlist = sorted(returnlist, key=operator.itemgetter('dist')) # sort stars according to dist

    return returnlist
