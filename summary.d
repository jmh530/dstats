/**Summary statistics such as mean, median, sum, variance, skewness, kurtosis.
 * Except for median and median absolute deviation, which cannot be calculated
 * online, all summary statistics have both an input range interface and an
 * output range interface.
 *
 * Bugs:  This whole module assumes that input will be reals or types implicitly
 *        convertible to real.  No allowances are made for user-defined numeric
 *        types such as BigInts.  This is necessary for simplicity.  However,
 *        if you have a function that converts your data to reals, most of
 *        these functions work with any input range, so you can simply map
 *        this function onto your range.
 *
 * Author:  David Simcha*/
 /*
 * License:
 * Boost Software License - Version 1.0 - August 17th, 2003
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */


module dstats.summary;

import std.algorithm, std.functional, std.conv, std.string, std.range,
       std.array;

import dstats.sort, dstats.base, dstats.alloc;

version(unittest) {
    import std.stdio, std.random, std.algorithm, std.conv;

    void main() {
    }
}

/**Finds median of an input range in O(N) time on average.  In the case of an
 * even number of elements, the mean of the two middle elements is returned.
 * This is a convenience founction designed specifically for numeric types,
 * where the averaging of the two middle elements is desired.  A more general
 * selection algorithm that can handle any type with a total ordering, as well
 * as selecting any position in the ordering, can be found at
 * dstats.sort.quickSelect() and dstats.sort.partitionK().
 * Allocates memory, does not reorder input data.*/
real median(T)(T data)
if(realInput!(T)) {
    // Allocate once on TempAlloc if possible, i.e. if we know the length.
    // This can be done on TempAlloc.  Otherwise, have to use GC heap
    // and appending.
    auto dataDup = tempdup(data);
    scope(exit) TempAlloc.free;
    return medianPartition(dataDup);
}

/**Median finding as in median(), but will partition input data such that
 * elements less than the median will have smaller indices than that of the
 * median, and elements larger than the median will have larger indices than
 * that of the median. Useful both for its partititioning and to avoid
 * memory allocations.  Requires a random access range with swappable
 * elements.*/
real medianPartition(T)(T data)
if(isRandomAccessRange!(T) &&
   is(ElementType!(T) : real) &&
   hasSwappableElements!(T) &&
   dstats.base.hasLength!(T))
{
    if(data.length == 0) {
        return real.nan;
    }
    // Upper half of median in even length case is just the smallest element
    // with an index larger than the lower median, after the array is
    // partially sorted.
    if(data.length == 1) {
        return data[0];
    } else if(data.length & 1) {  //Is odd.
        return cast(real) partitionK(data, data.length / 2);
    } else {
        auto lower = partitionK(data, data.length / 2 - 1);
        auto upper = ElementType!(T).max;

        // Avoid requiring slicing to be supported.
        foreach(i; data.length / 2..data.length) {
            if(data[i] < upper) {
                upper = data[i];
            }
        }
        return lower * 0.5L + upper * 0.5L;
    }
}

unittest {
    float brainDeadMedian(float[] foo) {
        qsort(foo);
        if(foo.length & 1)
            return foo[$ / 2];
        return (foo[$ / 2] + foo[$ / 2 - 1]) / 2;
    }

    float[] test = new float[1000];
    uint upperBound, lowerBound;
    foreach(testNum; 0..1000) {
        foreach(ref e; test) {
            e = uniform(0f, 1000f);
        }
        do {
            upperBound = uniform(0u, test.length);
            lowerBound = uniform(0u, test.length);
        } while(lowerBound == upperBound);
        if(lowerBound > upperBound) {
            swap(lowerBound, upperBound);
        }
        auto quickRes = median(test[lowerBound..upperBound]);
        auto accurateRes = brainDeadMedian(test[lowerBound..upperBound]);

        // Off by some tiny fraction in even N case because of division.
        // No idea why, but it's too small a rounding error to care about.
        assert(approxEqual(quickRes, accurateRes));
    }

    // Make sure everything works with lowest common denominator range type.
    struct Count {
        uint num;
        uint upTo;
        uint front() {
            return num;
        }
        void popFront() {
            num++;
        }
        bool empty() {
            return num >= upTo;
        }
    }

    Count a;
    a.upTo = 100;
    assert(approxEqual(median(a), 49.5));
    writeln("Passed median unittest.");
}

/**Calculates the median absolute deviation of a dataset.  This is the median
 * of all absolute differences from the median of the dataset.
 *
 * Notes:  No bias correction is used in this implementation, since using
 * one would require assumptions about the underlying distribution of the data.
 */
real medianAbsDev(T)(T data)
if(realInput!(T)) {
    auto dataDup = tempdup(data);
    immutable med = medianPartition(dataDup);
    immutable len = dataDup.length;
    TempAlloc.free;

    real[] devs = newStack!real(len);

    size_t i = 0;
    foreach(elem; data) {
        devs[i++] = abs(med - elem);
    }
    auto ret = medianPartition(devs);
    TempAlloc.free;
    return ret;
}

unittest {
    assert(approxEqual(medianAbsDev([7,1,8,2,8,1,9,2,8,4,5,9].dup), 2.5L));
    assert(approxEqual(medianAbsDev([8,6,7,5,3,0,999].dup), 2.0L));
    writeln("Passed medianAbsDev unittest.");
}

/**Output range to calculate the mean online.  Getter for mean costs a branch to
 * check for N == 0.  This struct uses O(1) space and does *NOT* store the
 * individual elements.
 *
 * Note:  This struct can implicitly convert to the value of the mean.
 *
 * Examples:
 * ---
 * Mean summ;
 * summ.put(1);
 * summ.put(2);
 * summ.put(3);
 * summ.put(4);
 * summ.put(5);
 * assert(summ.mean == 3);
 * ---*/
struct Mean {
private:
    real result = 0;
    real k = 0;

public:
    /// Allow implicit casting to real, by returning the current mean.
    alias mean this;

    ///
    void put(real element) nothrow {
        result += (element - result) / ++k;
    }

    /**Adds the contents of rhs to this instance.
     *
     * Examples:
     * ---
     * Mean mean1, mean2, combined;
     * foreach(i; 0..5) {
     *     mean1.put(i);
     * }
     *
     * foreach(i; 5..10) {
     *     mean2.put(i);
     * }
     *
     * mean1.put(mean2);
     *
     * foreach(i; 0..10) {
     *     combined.put(i);
     * }
     *
     * assert(approxEqual(combined.mean, mean1.mean));
     * ---
     */
     void put(const ref typeof(this) rhs) nothrow {
         immutable totalN = k + rhs.k;
         result = result * (k / totalN) + rhs.result * (rhs,k / totalN);
         k = totalN;
     }

    ///
    real sum() const pure nothrow {
        return result * k;
    }

    ///
    real mean() const pure nothrow {
        return (k == 0) ? real.nan : result;
    }

    ///
    real N() const pure nothrow {
        return k;
    }

    ///
    string toString() const {
        return to!(string)(mean);
    }
}

/**Finds the arithmetic mean of any input range whose elements are implicitly
 * convertible to real.*/
Mean mean(T)(T data)
if(realIterable!(T)) {
    Mean meanCalc;
    foreach(element; data) {
        meanCalc.put(element);
    }
    return meanCalc;
}

///
struct GeometricMean {
private:
    Mean m;
public:
    ///Allow implicit casting to real, by returning current geometric mean.
    alias geoMean this;

    ///
    void put(real element) nothrow {
        m.put(log2(element));
    }

    /// Combine two GeometricMean's.
    void put(const ref typeof(this) rhs) nothrow {
        m.put(rhs.m);
    }

    ///
    real geoMean() const pure nothrow {
        return exp2(m.mean);
    }

    ///
    real N() const pure nothrow {
        return m.k;
    }

    ///
    string toString() const {
        return to!(string)(geoMean);
    }
}

///
real geometricMean(T)(T data)
if(realIterable!(T)) {
    GeometricMean m;
    foreach(elem; data) {
        m.put(elem);
    }
    return m.geoMean;
}

unittest {
    string[] data = ["1", "2", "3", "4", "5"];
    auto foo = map!(to!(uint, string))(data);

    auto result = geometricMean(map!(to!(uint, string))(data));
    assert(approxEqual(result, 2.60517));

    Mean mean1, mean2, combined;
    foreach(i; 0..5) {
      mean1.put(i);
    }

    foreach(i; 5..10) {
      mean2.put(i);
    }

    mean1.put(mean2);

    foreach(i; 0..10) {
      combined.put(i);
    }

    assert(approxEqual(combined.mean, mean1.mean));
    assert(combined.N == mean1.N);

    writeln("Passed geometricMean unittest.");
}


/**Finds the sum of an input range whose elements implicitly convert to real.
 * User has option of making U a different type than T to prevent overflows
 * on large array summing operations.  However, by default, return type is
 * T (same as input type).*/
U sum(T, U = Unqual!(IterType!(T)))(T data)
if(realIterable!(T)) {
    U sum = 0;
    foreach(value; data) {
        sum += value;
    }
    return sum;
}

unittest {
    assert(sum(cast(int[]) [1,2,3,4,5])==15);
    assert(approxEqual( sum(cast(int[]) [40.0, 40.1, 5.2]), 85.3));
    assert(mean(cast(int[]) [1,2,3]) == 2);
    assert(mean(cast(int[]) [1.0, 2.0, 3.0]) == 2.0);
    assert(mean([1, 2, 5, 10, 17][]) == 7);
    assert(mean([1, 2, 5, 10, 17][]).sum == 35);
    writeln("Passed sum/mean unittest.");
}


/**Output range to compute mean, stdev, variance online.  Getter methods
 * for stdev, var cost a few floating point ops.  Getter for mean costs
 * a single branch to check for N == 0.  Relatively expensive floating point
 * ops, if you only need mean, try Mean.  This struct uses O(1) space and
 * does *NOT* store the individual elements.
 *
 * Note:  This struct can implicitly convert to a Mean struct.
 *
 * Examples:
 * ---
 * MeanSD summ;
 * summ.put(1);
 * summ.put(2);
 * summ.put(3);
 * summ.put(4);
 * summ.put(5);
 * assert(summ.mean == 3);
 * assert(summ.stdev == sqrt(2.5));
 * assert(summ.var == 2.5);
 * ---*/
struct MeanSD {
private:
    real _mean = 0;
    real _var = 0;
    real _k = 0;
public:
    ///
    void put(real element) nothrow {
        real kNeg1 = 1.0L / ++_k;
        _var += (element * element - _var) * kNeg1;
        _mean += (element - _mean) * kNeg1;
    }

    /// Combine two MeanSD's.
    void put(const ref typeof(this) rhs) nothrow {
        immutable totalN = _k + rhs._k;
        _mean = _mean * (_k / totalN) + rhs._mean * (rhs._k / totalN);
        _var = _var * (_k / totalN) + rhs._var * (rhs._k / totalN);
        _k = totalN;
    }

    ///
    real sum() const pure nothrow {
        return _k * _mean;
    }

    ///
    real mean() const pure nothrow {
        return (_k == 0) ? real.nan : _mean;
    }

    ///
    real stdev() const pure nothrow {
        return sqrt(var);
    }

    ///
    real var() const pure nothrow {
        return (_k < 2) ? real.nan : (_var - _mean * _mean) * (_k / (_k - 1));
    }

    // Undocumented on purpose b/c it's for internal use only.
    real mse() const pure nothrow {
        return (_k < 2) ? real.nan : (_var - _mean * _mean);
    }

    ///
    real N() const pure nothrow {
        return _k;
    }

    /**Converts this struct to a Mean struct.  Also called when an
     * implicit conversion via alias this takes place.
     */
    Mean toMean() const pure nothrow {
        return Mean(_mean, _k);
    }

    ///
    string toString() const {
        return text("N = ", cast(ulong) _k, "\nMean = ", mean, "\nVariance = ",
               var, "\nStdev = ", stdev);
    }
}

/**Convenience function that puts all elements of data into a MeanSD struct,
 * then returns this struct.*/
MeanSD meanStdev(T)(T data)
if(realIterable!(T)) {
    MeanSD ret;
    foreach(elem; data) {
        ret.put(elem);
    }
    return ret;
}

/**Finds the variance of an input range with members implicitly convertible
 * to reals.*/
real variance(T)(T data)
if(realIterable!(T)) {
    return meanStdev(data).var;
}

/**Calculate the standard deviation of an input range with members
 * implicitly converitble to real.*/
real stdev(T)(T data)
if(realIterable!(T)) {
    return meanStdev(data).stdev;
}

unittest {
    auto res = meanStdev(cast(int[]) [3, 1, 4, 5]);
    assert(approxEqual(res.stdev, 1.7078));
    assert(approxEqual(res.mean, 3.25));
    res = meanStdev(cast(double[]) [1.0, 2.0, 3.0, 4.0, 5.0]);
    assert(approxEqual(res.stdev, 1.5811));
    assert(approxEqual(res.mean, 3));
    assert(approxEqual(res.sum, 15));

    MeanSD mean1, mean2, combined;
    foreach(i; 0..5) {
      mean1.put(i);
    }

    foreach(i; 5..10) {
      mean2.put(i);
    }

    mean1.put(mean2);

    foreach(i; 0..10) {
      combined.put(i);
    }

    assert(approxEqual(combined.mean, mean1.mean));
    assert(approxEqual(combined.stdev, mean1.stdev));
    assert(combined.N == mean1.N);

    writefln("Passed variance/standard deviation unittest.");
}

/**Output range to compute mean, stdev, variance, skewness, kurtosis, min, and
 * max online. Using this struct is relatively expensive, so if you just need
 * mean and/or stdev, try MeanSD or Mean. Getter methods for stdev,
 * var cost a few floating point ops.  Getter for mean costs a single branch to
 * check for N == 0.  Getters for skewness and kurtosis cost a whole bunch of
 * floating point ops.  This struct uses O(1) space and does *NOT* store the
 * individual elements.
 *
 * Note:  This struct can implicitly convert to a MeanSD.
 *
 * Examples:
 * ---
 * Summary summ;
 * summ.put(1);
 * summ.put(2);
 * summ.put(3);
 * summ.put(4);
 * summ.put(5);
 * assert(summ.N == 5);
 * assert(summ.mean == 3);
 * assert(summ.stdev == sqrt(2.5));
 * assert(summ.var == 2.5);
 * assert(approxEqual(summ.kurtosis, -1.9120));
 * assert(summ.min == 1);
 * assert(summ.max == 5);
 * assert(summ.sum == 15);
 * ---*/
struct Summary {
private:
    real _mean = 0;
    real _m2 = 0;
    real _m3 = 0;
    real _m4 = 0;
    real _k = 0;
    real _min = real.infinity;
    real _max = -real.infinity;
public:
    ///
    void put(real element) nothrow {
        immutable real kNeg1 = 1.0L / ++_k;
        _min = (element < _min) ? element : _min;
        _max = (element > _max) ? element : _max;
        _mean += (element - _mean) * kNeg1;
        _m2 += (element * element - _m2) * kNeg1;
        _m3 += (element * element * element - _m3) * kNeg1;
        _m4 += (element * element * element * element - _m4) * kNeg1;
    }

    /// Combine two Summary's.
    void put(const ref typeof(this) rhs) nothrow {
        immutable totalN = _k + rhs._k;
        _mean = _mean * (_k / totalN) + rhs._mean * (rhs._k / totalN);
        _m2 = _m2 * (_k / totalN) + rhs._m2 * (rhs._k / totalN);
        _m3 = _m3 * (_k / totalN) + rhs._m3 * (rhs._k / totalN);
        _m4 = _m4 * (_k / totalN) + rhs._m4 * (rhs._k / totalN);
        _min = (_min < rhs._min) ? _min : rhs._min;
        _max = (_max > rhs._max) ? _max : rhs._max;
        _k = totalN;
    }

    ///
    real sum() const pure nothrow {
        return _mean * _k;
    }

    ///
    real mean() const pure nothrow {
        return (_k == 0) ? real.nan : _mean;
    }

    ///
    real stdev() const pure nothrow {
        return sqrt(var);
    }

    ///
    real var() const pure nothrow {
        return (_k == 0) ? real.nan : (_m2 - _mean * _mean) * (_k / (_k - 1));
    }

    ///
    real skewness() const pure nothrow {
        real var = _m2 - _mean * _mean;
        real numerator = _m3 - 3 * _mean * _m2 + 2 * _mean * _mean * _mean;

        // Raising var to the power of 1.5.  Non-obvious method is faater than
        // calling pow and allows this funciton to be pure nothrow.
        real sd = sqrt(var);
        real var15 = sd * sd * sd;
        return numerator / var15;
    }

    ///
    real kurtosis() const pure nothrow {
        real mean4 = mean * mean;
        mean4 *= mean4;
        real vari = _m2 - _mean * _mean;
        return (_m4 - 4 * _mean * _m3 + 6 * _mean * _mean * _m2 - 3 * mean4) /
               (vari * vari) - 3;
    }

    ///
    real N() const pure nothrow {
        return _k;
    }

    ///
    real min() const pure nothrow {
        return _min;
    }

    ///
    real max() const pure nothrow {
        return _max;
    }

    /**Converts this struct to a MeanSD.  Called via alias this when an
     * implicit conversion is attetmpted.
     */
    MeanSD toMeanSD() const pure nothrow {
        return MeanSD(_mean, _m2, _k);
    }

    alias toMeanSD this;

    ///
    string toString() const {
        return text("N = ", roundTo!long(_k),
                  "\nMean = ", mean,
                  "\nVariance = ", var,
                  "\nStdev = ", stdev,
                  "\nSkewness = ", skewness,
                  "\nKurtosis = ", kurtosis,
                  "\nMin = ", _min,
                  "\nMax = ", _max);
    }
}

unittest {
    // Everything else is tested indirectly through kurtosis, skewness.  Test
    // put(typeof(this)).

    Summary mean1, mean2, combined;
    foreach(i; 0..5) {
      mean1.put(i);
    }

    foreach(i; 5..10) {
      mean2.put(i);
    }

    mean1.put(mean2);

    foreach(i; 0..10) {
      combined.put(i);
    }

    foreach(ti, elem; mean1.tupleof) {
        assert(approxEqual(elem, combined.tupleof[ti]));
    }
}

/**Excess kurtosis relative to normal distribution.  High kurtosis means that
 * the variance is due to infrequent, large deviations from the mean.  Low
 * kurtosis means that the variance is due to frequent, small deviations from
 * the mean.  The normal distribution is defined as having kurtosis of 0.
 * Input must be an input range with elements implicitly convertible to real.*/
real kurtosis(T)(T data)
if(realIterable!(T)) {
    Summary kCalc;
    foreach(elem; data) {
        kCalc.put(elem);
    }
    return kCalc.kurtosis;
}

unittest {
    // Values from Matlab.
    assert(approxEqual(kurtosis([1, 1, 1, 1, 10].dup), 0.25));
    assert(approxEqual(kurtosis([2.5, 3.5, 4.5, 5.5].dup), -1.36));
    assert(approxEqual(kurtosis([1,2,2,2,2,2,100].dup), 2.1657));
    writefln("Passed kurtosis unittest.");
}

/**Skewness is a measure of symmetry of a distribution.  Positive skewness
 * means that the right tail is longer/fatter than the left tail.  Negative
 * skewness means the left tail is longer/fatter than the right tail.  Zero
 * skewness indicates a symmetrical distribution.  Input must be an input
 * range with elements implicitly convertible to real.*/
real skewness(T)(T data)
if(realIterable!(T)) {
    Summary sCalc;
    foreach(elem; data) {
        sCalc.put(elem);
    }
    return sCalc.skewness;
}

unittest {
    // Values from Octave.
    assert(approxEqual(skewness([1,2,3,4,5].dup), 0));
    assert(approxEqual(skewness([3,1,4,1,5,9,2,6,5].dup), 0.5443));
    assert(approxEqual(skewness([2,7,1,8,2,8,1,8,2,8,4,5,9].dup), -0.0866));

    // Test handling of ranges that are not arrays.
    string[] stringy = ["3", "1", "4", "1", "5", "9", "2", "6", "5"];
    auto intified = map!(to!(int, string))(stringy);
    assert(approxEqual(skewness(intified), 0.5443));
    writeln("Passed skewness test.");
}

/**Convenience function.  Puts all elements of data into a Summary struct,
 * and returns this struct.*/
Summary summary(T)(T data)
if(realIterable!(T)) {
    Summary summ;
    foreach(elem; data) {
        summ.put(elem);
    }
    return summ;
}
// Just a convenience function for a well-tested struct.  No unittest really
// necessary.  (Famous last words.)

///
struct ZScore(T) if(isForwardRange!(T) && is(ElementType!(T) : real)) {
private:
    T range;
    real mean;
    real sdNeg1;

    real z(real elem) {
        return (elem - mean) * sdNeg1;
    }

public:
    this(T range) {
        this.range = range;
        auto msd = meanStdev(range);
        this.mean = msd.mean;
        this.sdNeg1 = 1.0L / msd.stdev;
    }

    this(T range, real mean, real sd) {
        this.range = range;
        this.mean = mean;
        this.sdNeg1 = 1.0L / sd;
    }

    ///
    real front() {
        return z(range.front);
    }

    ///
    void popFront() {
        range.popFront;
    }

    ///
    bool empty() {
        return range.empty;
    }

    static if(isRandomAccessRange!(T)) {
        ///
        real opIndex(size_t index) {
            return z(range[index]);
        }
    }

    static if(isBidirectionalRange!(T)) {
        ///
        real back() {
            return z(range.back);
        }

        ///
        void popBack() {
            range.popBack;
        }
    }

    static if(dstats.base.hasLength!(T)) {
        ///
        size_t length() {
            return range.length;
        }
    }
}

/**Returns a range with whatever properties T has (forward range, random
 * access range, bidirectional range, hasLength, etc.),
 * of the z-scores of the underlying
 * range.  A z-score of an element in a range is defined as
 * (element - mean(range)) / stdev(range).
 *
 * Notes:
 *
 * If the data contained in the range is a sample of a larger population,
 * rather than an entire population, then technically, the results output
 * from the ZScore range are T statistics, not Z statistics.  This is because
 * the sample mean and standard deviation are only estimates of the population
 * parameters.  This does not affect the mechanics of using this range,
 * but it does affect the interpretation of its output.
 *
 * Accessing elements of this range is fairly expensive, as a
 * floating point multiply is involved.  Also, constructing this range is
 * costly, as the entire input range has to be iterated over to find the
 * mean and standard deviation.
 */
ZScore!(T) zScore(T)(T range)
if(isForwardRange!(T) && realInput!(T)) {
    return ZScore!(T)(range);
}

/**Allows the construction of a ZScore range with precomputed mean and
 * stdev.
 */
ZScore!(T) zScore(T)(T range, real mean, real sd)
if(isForwardRange!(T) && realInput!(T)) {
    return ZScore!(T)(range, mean, sd);
}

unittest {
    int[] arr = [1,2,3,4,5];
    auto m = mean(arr);
    auto sd = stdev(arr);
    auto z = zScore(arr);

    size_t pos = 0;
    foreach(elem; z) {
        assert(elem == (arr[pos++] - m) / sd);
    }

    assert(z.length == 5);
    foreach(i; 0..z.length) {
        assert(z[i] == (arr[i] - m) / sd);
    }
    writeln("Passed zScore test.");
}



// Verify that there are no TempAlloc memory leaks anywhere in the code covered
// by the unittest.  This should always be the last unittest of the module.
unittest {
    auto TAState = TempAlloc.getState;
    assert(TAState.used == 0);
    assert(TAState.nblocks < 2);
}
