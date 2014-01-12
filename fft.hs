-- FFT and Spectrum Analyzer
-- Bohumir Zamecnik <bohumir [at] zamecnik [dot] org>

import Complex
import Array

--import Debug.Trace

-- References:
-- FFT algorithms used here are implemented according the book
-- Cormen, Leiserson, Rivest, Stein: Introduction to Algorithms, 2nd Ed., 2001

-- FFT Algorithms ----------------------------------------------------

-- FFT recursively -----------------------------------------

-- FFT recursively on complex lists --------------

-- fft & ifft are just wrappers of fft'
-- fft :: (RealFloat a) => [Complex a] -> [Complex a]
fft []  = []
fft xs =  fft' xs' (imagroots n')
	where
		xs' = padRightToPowerOf2 xs
		n' = length xs'
		n = length xs

-- Inverse FFT recursively on complex lists ------

-- fft :: (RealFloat a) => [Complex a] -> [Complex a]
ifft []  = []
ifft xs = map ((/(fromIntegral n)).conjugate) (fft' xs' (imagroots n))
	where
		xs' = map conjugate (padRightToPowerOf2 xs)
		n = length xs'

-- FFT recursively on lists - core ---------------

-- fft' :: Num a => [a] -> [a] -> [a]
fft' [] _  = []
fft' [x] _ = [x]
--fft' _ [] = []
fft' xs roots = (zipWith (+) evenxs $ products)
                 ++ (zipWith (-) evenxs $ products)
	where
		oddxs = fft' (odds xs) evenroots
		evenxs = fft' (evens xs) evenroots
		evenroots = evens roots
		products = zipWith (*) oddxs roots

-- FFT recursively on real lists -----------------

-- fftR :: RealFloat a => [a] -> [Complex a]
fftR xs = fftR' fft

-- fftR :: RealFloat a => [a] -> [Complex a]
ifftR xs = fftR' ifft

-- fftR' :: RealFloat a => ([Complex a] -> b) -> [a] -> b
fftR' f xs = f (map (:+0) xs)


-- FFT iteratively -----------------------------------------

-- FFT iteratively on complex arrays -------------

-- Remark: FOR cycles are implemented as recursive calls
-- TODO: Make array handling more memory efficient!

-- fftiter :: (RealFloat a, Ix b) => Array b (Complex a)
--             -> Array Int (Complex a)
fftiterArray a = fftiter' 1 (log2 (fromIntegral n)) roots a'
	where
		a' = listArray (0,n-1) (bitReverseList padded)
		padded = padRightToPowerOf2 (elems a)
		n = length padded
		--n = arrayLength a' -- another implementation
		roots = imagrootsArray n

-- fftiter' :: (Integral a, Num b, Ix a, Integral c) => c -> c -> Array a b
--              -> Array a b -> Array a b
fftiter' s to roots a
	| s <= to =
--		trace ("fftiter' "--a="++(show (elems a))
--			++" s="++(show s)
--			++" to="++(show to))
		(fftiter' (s+1) to roots (fftiter'' 0 s m roots a))
--	| otherwise = trace "fftiter' return" a
	| otherwise = a
	where
		m = 2^s

-- fftiter'' :: (Ix a, Num b, Integral a) => a -> c -> a -> Array a b
--               -> Array a b -> Array a b
fftiter'' k s m roots a
	| k < n =
--		trace ("fftiter''' "--a="++(show (elems a))
--			++" k="++(show k)
--			++" s="++(show s))
		(fftiter'' (k+m) s m roots (fftiter''' 0 s k m roots a))
--	| otherwise = trace "fftiter' return" a
	| otherwise = a
	where
		n = arrayLength a

-- fftiter''' :: (Integral a, Num b, Ix a) => a -> c -> a -> a -> Array a b
--                -> Array a b -> Array a b
fftiter''' j s k m roots a
	| j < (m `div` 2) =
--		trace ("fftiter'' "--a="++(show (elems a))
--			++" j="++(show j)
--			++" m="++(show m))
		(fftiter''' (j+1) s k m roots (butterflyArray k j m roots a))
--	| otherwise = trace "fftiter' return" a
	| otherwise = a
	where
		n = arrayLength a

-- butterfly :: (Ix a, Num b, Integral a) => a -> a -> a -> Array a b
--               -> Array a b -> Array a b
butterflyArray k j m roots a
	| (k+j+m2) < n =
--		trace ("butterfly " ++" k="++(show k) ++" m="++(show m))
		(a // [(k+j,u+t),(k+j+m2,u-t)])
	| otherwise = a
	where
		t = roots!(j*q) * a!(k+j+m2)
		u = a!(k+j)
		q = n `div` m
		n = arrayLength a
		m2 = m `div` 2

-- Inverse FFT iteratively on complex arrays -----

-- ifftiter :: (Ix a, RealFloat b) => Array a (Complex b)
--              -> Array Int (Complex b)
--ifftiter a = listArray (0,n-1) outlist
--	where
--		inlist = map conjugate (padRightToPowerOf2 (elems a))
--		n = length inlist
--		outlist = map ((/(fromIntegral n)).conjugate) (fftiterList inlist)

-- FFT iteratively on lists ---------------------------

fftiter = fftiterList
ifftiter = ifftiterList

-- fftiterList :: RealFloat a => [Complex a] -> [Complex a]
fftiterList xs = unwrap (fftiterList' 2 n xs')
	where
		xs' = wrap (bitReverseList padded)
		padded = padRightToPowerOf2 xs
		n = length padded

-- fftiterList' :: RealFloat a => Int -> Int -> [[Complex a]] -> [[Complex a]]
fftiterList' s to xs
	| s <= to = fftiterList' (s*2) to (bfmap (imagroots s) xs)
	| otherwise = xs

-- ifftiterList :: RealFloat a => [Complex a] -> [Complex a]
ifftiterList xs = map ((/(fromIntegral n)).conjugate) (fftiterList xs')
	where
		xs' = map conjugate (padRightToPowerOf2 xs)
		n = length xs'

-- helper functions for fftiterList

-- bfmap :: Num a => [a] -> [[a]] -> [[a]]
bfmap _ [] = []
bfmap omegas (x1:x2:xs) = (butterflyZip omegas x1 x2 [] []):(bfmap omegas xs)

--bfmap omegas xs = butterflyZip omegas (take n2 xs) (drop n2 xs)
--	where
--		n2 = (length xs) `div` 2

-- butterflyZip :: Num a => [a] -> [a] -> [a] -> [a] -> [a] -> [a]
butterflyZip _ [] [] a1 a2 = concat [reverse a1,reverse a2]
butterflyZip (o:omegas) (x:xs) (y:ys) a1 a2
	= butterflyZip omegas xs ys ((fst b):a1) ((snd b):a2)
	where
		b = butterfly o x y

-- butterfly operation
-- butterfly :: Num a => a -> a -> a -> (a,a)
butterfly omega a b = (a+t,a-t)
	where
		t = omega*b

-- wrap :: [a] -> [[a]]
wrap = map (\x->[x])
-- unwrap :: [a] -> a
unwrap [x] = x

-- FFT iteratively on complex lists --------------

-- fftiterList :: [Complex Double] -> [Complex Double]
--fftiterList = fftiterList' (fftiter)

-- ifftiterList :: [Complex Double] -> [Complex Double]
--ifftiterList = fftiterList' (ifftiter)

-- fftiterList' :: RealFloat a => [Complex a] -> [Complex a]
--fftiterList' f xs = elems (f (listArray (0,(length xs)-1) xs))

-- FFT iteratively on real lists -----------------

-- fftiterListR :: [Double] -> [Complex Double]
--fftiterListR = fftiterListR' (fftiterList)

-- ifftiterListR :: [Double] -> [Complex Double]
--ifftiterListR = fftiterListR' (ifftiterList)

-- fftiterListR' :: RealFloat a => [a] -> [Complex a]
--fftiterListR' f xs = f (map (:+0) xs)


-- DFT naively ---------------------------------------------

-- DFT on complex lists --------------------------

-- dft :: RealFloat a => [Complex a] -> [Complex a]
dft [x] = [x]
dft xs = dftnaive' xs' (imagroots n)
	where
		xs' = padRightToPowerOf2 xs
		n = length xs'

-- Inverse DFT on complex lists ------------------

-- idft :: RealFloat a => [Complex a] -> [Complex a]
idft xs = map ((/(fromIntegral n)).conjugate) (dftnaive' xs' (imagroots n))
	where
	xs' = map conjugate (padRightToPowerOf2 xs)
	n = length xs'

-- DFT - core ------------------------------------

-- dftnaive' :: Num a => [a] -> [a] -> [a]
dftnaive' [x] _ = [x]
dftnaive' xs roots = [ sum [x * omega^k|(x,omega)<-zip xs roots] |k<-[0..n-1]]
	where
		n = length xs


-- Helper functions for FFT etc. ---------------------------

-- Base 2 logarithm ------------------------------

-- log2 :: (RealFrac a, Integral b, Floating a) => a -> b
log2 n = ceiling (logBase 2 n)
--log2 n = ceiling (((log n)/(log 2))) -- another implementation

-- Next power of 2 <= n --------------------------

-- nextPowerOf2 :: (Num a, Integral b) => b -> a
nextPowerOf2 0 = 1
nextPowerOf2 n = 2 ^ (log2 (fromIntegral n))

-- Is power of 2? --------------------------------

-- isPowerOf2 :: Integral a => a -> Bool
isPowerOf2 n = n == nextPowerOf2 n
-- another implementation:
--isPowerOf2 1 = True
--isPowerOf2 2 = True
--isPowerOf2 n = ((n `mod` 2) == 0) && (isPowerOf2 (n `div` 2))

-- Length of array (eg. number of its elements) --

-- arrayLength :: (Num a, Ix a) => Array a b -> a
arrayLength a = 1 + (snd (bounds a)) - (fst (bounds a))

-- Pad list of length n with zeros to the nearest 2^k >= n --

-- padRightToPowerOf2 :: Num a => [a] -> [a]
padRightToPowerOf2 [] = []
padRightToPowerOf2 xs = xs ++ (replicate ((nextPowerOf2 n) - n) 0)
	where n = length xs

-- Select even items from the list ---------------

-- evens :: [a] -> [a]
evens xs = [x | (i,x) <- zip [0..] xs, even i]

-- Select odd items from the list ----------------

-- odds :: [a] -> [a]
odds xs = [x | (i,x) <- zip [0..] xs, odd i]

-- Sort a list in Bit Reverse order using recursion --
--	eg. the bits in indices are reversed 

-- bitReverseList :: [a] -> [a]
bitReverseList :: [a] -> [a]
bitReverseList [] = []
bitReverseList [x] = [x]
bitReverseList xs = bitReverseList (evens xs) ++ bitReverseList (odds xs)

-- Sort an array in Bit Reverse order ------------

-- bitReverseArray :: Ix a => Array a b -> Array a b
bitReverseArray a = listArray (bounds a) (bitReverseList (elems a))

-- Reverse append (for DCT) ----------------------
--		Trim first and last item, reverse it
--		and append to the original list.
--	eg.:	[a,b,c,d] -> [a,b,c,d,c,b]

-- revappend :: [a] -> [a]
revappend [x]  = [x]
revappend (x:xs) = x:xs ++ (reverse (init xs))

--	another implementation:
--revappend []  = []
--revappend (x:xs) = revappend' xs (x:xs)
--	where
--		revappend' [_] acc = acc
--		revappend' (x:xs) acc = (revappend' xs acc)++[x]


-- Complex roots of unity, eg. n-th roots of 1 --------

-- List of n-th roots of imaginary unit ----------

-- imagroots :: RealFloat a => Int -> [Complex a]
imagroots n = imagroots' n (-1)

-- List of inverse n-th roots of imaginary unit --

-- imagrootsinv :: RealFloat a => Int -> [Complex a]
imagrootsinv n = imagroots' n 1

-- Array of n-th roots of the complex unit -------

-- imagrootsArray :: RealFloat a => Int -> Array Int (Complex a)
imagrootsArray n = listArray (0,n-1) (imagroots' n (-1))

-- Array of inverse n-th roots of complex unit ---

-- imagrootsinvArray :: RealFloat a => Int -> Array Int (Complex a)
imagrootsinvArray n = listArray (0,n-1) (imagroots' n 1)

--	another implementation:
--imagrootsinv n = (head (roots)):(reverse (tail (roots)))
--	where
--		roots = imagroots n

-- The core for computing n-th roots of imaginary unit --

-- imagroots' :: RealFloat a => Int -> a -> [Complex a]
imagroots' n sign = take n (iterate (*omega) 1)
	where
		-- omega is the basic n-th root
		omega = cis (sign*2*pi/fromIntegral n)


-- Epsilon rounding ------------------------------
--	if |x| < epsilon, treat it as zero

-- epsilon :: RealFrac a => a -> a
epsilon x | abs(roundx - x) < eps = roundx
	| otherwise = x
	where
	eps = 1e-9
	roundx = fromIntegral (round x)

-- The same as epsilon, but for both components of the complex number --
-- Remark: very useful!

-- epsiloncplx :: RealFloat a => Complex a -> Complex a
epsiloncplx (re:+im) = epsilon(re):+epsilon(im)

rnd = map epsiloncplx

-- Show a complex number (unused)
--showcplx (re:+im) = show re "+" im "i"

-- Convert a list of real numbers to a list of complex numbers --

-- cplxList :: [Double] -> [Complex Double]
cplxList = map (:+0)


-- Spectrum Analyzer -------------------------------------------------

-- Currently plots spectrum graphs with following features:
--	* in textmode (using bars of '*' characters)
--  * power scale is logarithmic (in dB-like units)
--	* frequency scale is linear

-- Calcutale real "power" from a complex number --
-- usage: map power [1:+2, ...]

-- power :: RealFloat a => Complex a -> a
power (r:+i) = r^2 + i^2

-- Calculate Y position in spectrogram from given power --
-- Remark: constants are taken, such that we have 100dB/80 character terminal
yheight = 80 -- points (pixels or chars)
ybels = 10.0 -- bels on scale

-- ploty :: Double -> Double
ploty 0 = 0
ploty p = ((logBase 10 p) * yscale)
	where
		yscale = (fromInteger yheight) / ybels

-- Plotting in text mode -----------------------------------

-- Make bar long x using '*' characters ----------

-- mkBar :: Integral a => a -> [Char]
mkBar x = (replicate (fromIntegral x) '*')++"\n"
-- alternative implementation:
--mkBar 0 = "\n"
--mkBar x = "*" ++ mkBr (x-1)

-- Draw bars on stdout given list of their lengths --

-- drawBars :: Integral a => [a] -> IO ()
drawBars xs = foldr (>>) (return ()) (map putStr (map mkBar xs))

-- Draw bars into a string (unused)
--mkBars xs = foldr (++) [] (map mkBar xs)

-- Calculate spectrum graph given a list ---------
-- with real values of input signal

-- spectrum :: [Double] -> [Double]
spectrum xs = take (ceiling (len*0.5)) s
	where
		s = map power (fftiterList (cplxList xs))
		len = fromIntegral (length s)

-- Plot spectrum graph using '*' chars -----------
-- Power scale is logarithmic (in dB-like units)

-- plotSpectrum :: [Double] -> IO ()
plotSpectrum xs = drawBars (map (ceiling.ploty) (spectrum xs))

-- Take input signal stream and use windowing method to calculate and plot
-- spectrum graph.
-- size: size of the window

-- processSignal :: Int -> [Double] -> IO ()
processSignal _ [] = return ()
processSignal size xs = (plotSpectrum (take size xs))
                         >> (processSignal size (drop size xs))

-- Tests -------------------------------------------------------------

-- FFT vs. DFT
-- Test if recursive FFT and naive DFT give the same result
testcompare fft1 fft2 list = (\l->foldr (&&) True (zipped l)) list
	where
		zipped l = (zipWith (==) (rnd (fft1 l)) (rnd (fft2 l)))

testfftdft = testcompare fft dft (imagroots 128)
testfftiter = testcompare fft fftiterList list4
testifftiter = testcompare id (ifftiter.fftiter) list1


-- some testing data
list1 = cplxList (replicate 1024 1)
list2 = cplxList [1..1024]
list3 = imagroots 1024
list4 = take 1024 (cycle (imagroots 512))

list5 = cplxList (replicate 256 1)
list6 = cplxList [1..256]
list7 = imagroots 256
list8 = take 256 (cycle (imagroots 128))

-- Test FFT
testfft1 = fft list1
testfft2 = fft list2
testfft3 = map epsiloncplx (fft list3)
testfft4 = map epsiloncplx (fft list4)

-- Test Inverse FFT
testifft1 = ifft (fft list1)
testifft2 = map epsiloncplx (ifft (fft list2))
testifft3 = map epsiloncplx (ifft (fft list3))
testifft4 = map epsiloncplx (ifft (fft list4))

-- FFT iteratively
-- Arrays
--testfftiterArray1 = testfftarray' fftiterArray list5
--testfftiterArray2 = testfftarray' fftiterArray list6
--testfftiterArray3 = testfftarray' fftiterArray list7
--testfftiterArray4 = testfftarray' fftiterArray list8

--testifftiter1 = testfftarray' (ifftiter.fftiter) list5
--testifftiter2 = testfftarray' (ifftiter.fftiter) list6
--testifftiter3 = testfftarray' (ifftiter.fftiter) list7
--testifftiter4 = testfftarray' (ifftiter.fftiter) list8

testfftarray' f xs = f (listArray (0,(length xs)-1) xs)

-- Lists
testfftiterList1 =  fftiterList list5
testfftiterList2 =  fftiterList list6
testfftiterList3 =  fftiterList list7
testfftiterList4 =  fftiterList list8

-- generate a sinus signal
singen freq len = map (sin) (map ((2*pi*freq/((fromIntegral len)))*) (take len [0..]))

testsingen = (\scale->map (ceiling.(scale+).(scale*)) (singen 2.0 32)) 16

---- test plotting
testpl1 = drawBars testsingen

---- test spectral analyzing
--testpl2 = drawBars (map (ceiling.ploty.power) (fft (cplxList (singen 1 32))))

-- plain sinus signal
testsp1 = plotSpectrum (singen 1 32)
-- two different sinusoids added
testsp2 = plotSpectrum (zipWith (+) (singen 1.0 32) (map (0.5*) (singen 5.0 32)))
-- sawtooth signal
testsp3 = plotSpectrum (take 31 (cycle ((take 4 [0..]) ++ (take 4 [4,3..]))))
-- square signal
testsp4 = plotSpectrum (take 31 (cycle ((replicate 4 0) ++ (replicate 4 1))))

-- text processing signal
testProcess = processSignal 32 [1..70]

