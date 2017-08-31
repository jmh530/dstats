/**Provides an alternate interface to access the dstats modules
 * dstats.distrib and dstats.random. This is not intended to be a replacement
 * of those modules, but rather an alternate approach to using those, rather it
 * offers an interface similar to Julia's distributions.jl.
 *
 * Each distribution has an associated struct that has the following member
 * functiosn defined, if they were available in the aforementioned dstats
 * modules:
 *      pdf or pmf - probability density/mass function
 *      density - either the probability density or probability mass function
 *      cdf - cumulative density function
 *      icdf - inverse cumulative density function
 *      cdfr - one minus the cumulative density function
 *      icdfr - inverse cdfr
 *      rand - random number generator
 *
 * For instance, the normal distribution is defined as the struct Normal
 * with the member functions: pdf, density, cdf, and icdf. This struct can be
 * used directly, as in Normal.cdf(0.5, 0, 1), which returns the cdf of a
 * standard normal evaluated at 0.5. The struct can also be used with
 * free-standing functions. For instance, one can call cdf!Normal(0.5, 0, 1) in
 * the same way as above.
 *
 * This module also provides implementations for dstats.random's randArray
 * and randRange to provide multiple random numbers.
 *
 * Finally, this module provides an alternative implementation of the
 * parametrize and paramFunctor functions in dstats.distrib, called
 * specificDistribution. The key difference is that the specificDistribution
 * function returns a functor (the instance of struct that can call any of the
 * above functions, assuming they exist for the underlying Distribution, with
 * specific parameters) that works for several functions, not just one. As the
 * behavior differs with respect to parametrize and paramFunctor, this
 * functionality is provided as a separate function.
 *
 * Authors:  John Hall
 */

module dstats.distributions;

import std.traits : isAggregateType;
import dstats.random : Random;

private template isMemberOf(T, string x)
    if (isAggregateType!T)
{
    import std.traits : hasMember;

    enum bool isMemberOf = hasMember!(T, x);
}

private void hasMemberCheck(alias T, string x)()
{
    static assert(isMemberOf!(T, x), T.stringof ~ " struct must have " ~ x ~
                                     " member to call " ~ x ~ " function");
}

private template hasRand(T)
    if (isAggregateType!T)
{
    enum bool hasRand = isMemberOf!(T, "rand");
}

private template hasCDF(T)
    if (isAggregateType!T)
{
    enum bool hasCDF = isMemberOf!(T, "cdf");
}

private template hasPDF(T)
    if (isAggregateType!T)
{
    enum bool hasPDF = isMemberOf!(T, "pdf");
}

private template hasPMF(T)
    if (isAggregateType!T)
{
    enum bool hasPMF = isMemberOf!(T, "pmf");
}

private template hasICDF(T)
    if (isAggregateType!T)
{
    enum bool hasICDF = isMemberOf!(T, "icdf");
}

private template hasDensity(T)
    if (isAggregateType!T)
{
    enum bool hasDensity = isMemberOf!(T, "density");
}

private template hasCDFR(T)
    if (isAggregateType!T)
{
    enum bool hasCDFR = isMemberOf!(T, "cdfr");
}

private template hasICDFR(T)
    if (isAggregateType!T)
{
    enum bool hasICDFR = isMemberOf!(T, "icdfr");
}

private template hasAnyDistrib(T)
    if (isAggregateType!T)
{
    enum bool hasAnyDistrib = hasCDF!T || hasPDF!T || hasPMF!T || hasICDF!T ||
                              hasDensity!T || hasCDFR!T || hasICDFR!T;
}

private template isTemplate(alias T)
{
    enum bool isTemplate = __traits(isTemplate, T);
}

private template hasAnyDistribTemplate(T)
    if (isAggregateType!T)
{
    static if (hasAnyDistrib!T)
    {
        static if (hasCDF!T && isTemplate!(T.cdf))
            enum bool hasAnyDistribTemplate = true;
        else static if (hasPDF!T && isTemplate!(T.pdf))
            enum bool hasAnyDistribTemplate = true;
        else static if (hasPMF!T && isTemplate!(T.pmf))
            enum bool hasAnyDistribTemplate = true;
        else static if (hasICDF!T && isTemplate!(T.icdf))
            enum bool hasAnyDistribTemplate = true;
        else static if (hasDensity!T && isTemplate!(T.density))
            enum bool hasAnyDistribTemplate = true;
        else static if (hasCDFR!T && isTemplate!(T.cdfr))
            enum bool hasAnyDistribTemplate = true;
        else static if (hasICDFR!T && isTemplate!(T.icdfr))
            enum bool hasAnyDistribTemplate = true;
        else
            enum bool hasAnyDistribTemplate = false;
    }
    else
    {
        enum bool hasAnyDistribTemplate = false;
    }
}

private template getProtection(string from, string member)
{
    mixin("static import " ~ from ~ ";");
    enum string getProtection =
                  mixin("__traits(getProtection, " ~ from ~ "." ~ member ~ ")");
}

@safe unittest
{
    assert(getProtection!("std.algorithm", "map") == "public");
}

// This currently also includes deprecated functions, so adjusted to just return
// true
private template isAvailable(string from, string member)
{
    enum bool isAvailable = true;
    /*
    enum string protect = getProtection!(from, member);
    static if (protect)
    {
        static if (protect == "public" ||
                   protect == "export")
        {
            enum bool isAvailable = true;
        }
        else
        {
            enum bool isAvailable = false;
        }
    }
    else
    {
        enum bool isAvailable = false;
    }
    */
}

private template getRName(string structName)
{
    static if (structName == "Normal")
        enum string getRName = "rNorm";
    else static if (structName == "LogNormal")
        enum string getRName = "rLogNorm";
    else static if (structName == "StudentsT")
        enum string getRName = "rStudentT";
    else
        enum string getRName = "r" ~ structName;
}

@safe unittest
{
    assert("rNorm" == getRName!"Normal");
    assert("rLogNorm" == getRName!"LogNormal");
    assert("rStudentT" == getRName!"StudentsT");
    assert("rBernoulli" == getRName!"Bernoulli");
}

private template getStructName(string rName)
{
    static if (rName == "rNorm")
        enum string getStructName = "Normal";
    else static if (rName == "rLogNorm")
        enum string getStructName = "LogNormal";
    else static if (rName == "rStudentT")
        enum string getStructName = "StudentsT";
    else
    {
        import std.algorithm.searching : findSplitAfter;
        enum string getStructName = findSplitAfter(rName, "r")[1];
    }
}

@safe unittest
{
    assert("Normal" == getStructName!"rNorm");
    assert("LogNormal" == getStructName!"rLogNorm");
    assert("StudentsT" == getStructName!"rStudentT");
    assert("Bernoulli" == getStructName!"rBernoulli");
}

private string genStructInternals(string funcName, string structName)()
{
    import dstats.distrib;
    import dstats.random;
    import std.array : appender;
    import std.string : toLower;
    import std.algorithm.searching : endsWith, startsWith, findSplitAfter;

    enum spaces = "    ";

    auto aliasBuf = appender!string();
    auto importDistBuf = appender!string();
    auto importRandBuf = appender!string();

    enum string invName = "inv" ~ structName;
    enum string rName = getRName!structName;

    bool hasDistrib = false;
    bool hasRand = false;

    importDistBuf.put(spaces);
    importDistBuf.put("import dstats.distrib : ");
    importRandBuf.put(spaces);
    importRandBuf.put("import dstats.random : ");

    foreach(member; __traits(allMembers, dstats.distrib))
    {
        static if (isAvailable!("dstats.distrib", member))
        {
            static if (startsWith(member, funcName))
            {
                enum string memberAfter = findSplitAfter(member, funcName)[1];
                enum string lowerMemberAfter = toLower(memberAfter);

                importDistBuf.put(member ~ ", ");

                aliasBuf.put(spaces);
                aliasBuf.put("alias " ~ lowerMemberAfter ~ " = "
                                                                ~ member ~ ";");
                aliasBuf.put("\n");

                if (hasDistrib == false)
                    hasDistrib = true;

                static if ((lowerMemberAfter == "pdf") ||
                           (lowerMemberAfter == "pmf"))
                {
                    aliasBuf.put(spaces);

                    aliasBuf.put("alias density = " ~ lowerMemberAfter ~ ";");
                    aliasBuf.put("\n");
                }
            }
            else static if (startsWith(member, invName))
            {
                enum string memberAfter = findSplitAfter(member, invName)[1];

                importDistBuf.put(member ~ ", ");

                aliasBuf.put(spaces);
                aliasBuf.put("alias i" ~ toLower(memberAfter) ~ " = "
                                                                ~ member ~ ";");
                aliasBuf.put("\n");
            }
        }
    }

    foreach(member; __traits(allMembers, dstats.random))
    {
        static if (isAvailable!("dstats.random", member))
        {
            static if (member == rName)
            {
                enum string memberAfter = findSplitAfter(member, rName)[1];

                importRandBuf.put(member ~ ";\n");

                aliasBuf.put(spaces);
                aliasBuf.put("alias rand = " ~ member ~ ";");
                aliasBuf.put("\n");

                if (hasRand == false)
                    hasRand = true;
            }
        }
    }

    if (hasRand == true || hasDistrib == true)
    {
        if (hasDistrib == true)
        {
            string preImportOut =
                           importDistBuf.data[0 .. ($ - (", ".length))] ~ ";\n";
            if (hasRand == true)
                return preImportOut ~ importRandBuf.data ~
                                        aliasBuf.data[0 .. ($ - ("\n").length)];
            else
                return preImportOut ~ aliasBuf.data[0 .. ($ - ("\n").length)];
        }
        else if (hasDistrib == false && hasRand == true)
        {
            return importRandBuf.data ~ aliasBuf.data[0 .. ($ - ("\n").length)];
        }
        else
        {
            assert(0, "Should not be here in generating internals for " ~
                                                                      funcName);
        }
    }
    else
    {
        assert(0,
            "No relevant functions in dstats.distrib or dstats.random for " ~
            funcName);
    }
}

@safe unittest
{
    enum string x = genStructInternals!("bernoulli", "Bernoulli");
    enum string y = "    import dstats.random : rBernoulli;\n" ~
                    "    alias rand = rBernoulli;";
    assert(x == y);
}

@safe unittest
{
    enum string x = genStructInternals!("dirichlet", "Dirichlet");
    enum string y = "    import dstats.distrib : dirichletPDF;\n" ~
                    "    alias pdf = dirichletPDF;\n" ~
                    "    alias density = pdf;";
    assert(x == y);
}

@safe unittest
{
    enum string x = genStructInternals!("normal", "Normal");
    enum string y = "    import dstats.distrib : normalPDF, normalCDF, " ~
                                               "normalCDFR, invNormalCDF;\n" ~
                    "    import dstats.random : rNorm;\n" ~
                    "    alias pdf = normalPDF;\n" ~
                    "    alias density = pdf;\n" ~
                    "    alias cdf = normalCDF;\n" ~
                    "    alias cdfr = normalCDFR;\n" ~
                    "    alias icdf = invNormalCDF;\n" ~
                    "    alias rand = rNorm;";
    assert(x == y);
}

private string toLowerFirst(string name)()
{
    import std.string : toLower;
    import std.conv : to;

    string firstLetter = name[0].toLower.to!string;
    return firstLetter ~ name[1 .. $];
}

@safe unittest
{
    enum string x = "LowerBlah";
    assert(toLowerFirst!x == "lowerBlah");
}

private string toUpperFirst(string name)()
{
    import std.string : toUpper;
    import std.conv : to;

    string firstLetter = name[0].toUpper.to!string;
    return firstLetter ~ name[1 .. $];
}

@safe unittest
{
    enum string x = "upperBlah";
    assert(toUpperFirst!x == "UpperBlah");
}

private template GenDistStruct(string name)
{
    const char[] GenDistStruct =
        "///" ~ "\n" ~
        "struct " ~ toUpperFirst!(name) ~ "\n" ~
        "{\n" ~
        genStructInternals!(name, toUpperFirst!(name)) ~ "\n" ~
        "}";
}

@safe unittest
{
    auto x = GenDistStruct!"bernoulli";
    assert(x == "///" ~ "\n" ~
            "struct Bernoulli" ~ "\n" ~
            "{\n" ~
            genStructInternals!("bernoulli", "Bernoulli") ~ "\n" ~
            "}");
}

private string GenDistStructs()
{
    import dstats.distrib;
    import dstats.random;
    import std.array : appender;
    import std.uni : isUpper;
    import std.algorithm.searching : startsWith, endsWith, canFind,
                                                findSplitBefore, findSplitAfter;

    string[__traits(allMembers, dstats.distrib).length +
           __traits(allMembers, dstats.random).length] createdStructs;
    size_t i;
    auto structsBuf = appender!string();

    foreach(member; __traits(allMembers, dstats.distrib))
    {
        static if (isAvailable!("dstats.distrib", member))
        {
            static if ((member.endsWith("PDF") ||
                        member.endsWith("PMF") ||
                        member.endsWith("CDF") ||
                        member.endsWith("CDFR")))
            {
                static if (member.endsWith("PDF"))
                    enum string memberBefore =
                                              findSplitBefore(member, "PDF")[0];
                else static if (member.endsWith("PMF"))
                    enum string memberBefore =
                                              findSplitBefore(member, "PMF")[0];
                else static if (member.endsWith("CDF"))
                    enum string memberBefore =
                                              findSplitBefore(member, "CDF")[0];
                else static if (member.endsWith("CDFR"))
                    enum string memberBefore =
                                             findSplitBefore(member, "CDFR")[0];

                static if (member.startsWith("inv"))
                    enum string newMember =
                          toLowerFirst!(findSplitAfter(memberBefore, "inv")[1]);
                else
                    enum string newMember = memberBefore;

                static if (member != "chiSqrCDF" &&
                           member != "chiSqrCDFR" &&
                           member != "invChiSqrCDFR" &&
                           member != "invChiSqCDFR")
                           //These are all deprecated, this is the only way I've
                           //figured out to fix it properly
                {
                    if (i == 0 ||
                               !(createdStructs[0 .. i].canFind(newMember)))
                    {
                        structsBuf.put(GenDistStruct!newMember);
                        structsBuf.put("\n");
                        createdStructs[i] = newMember;
                        i++;
                    }
                }
            }
        }
    }

    foreach(member; __traits(allMembers, dstats.random))
    {
        static if (isAvailable!("dstats.random", member))
        {
            static if (member.startsWith("r") && isUpper(member[1]))
            {
                enum string rName = getStructName!(member);
                enum string newMember = toLowerFirst!(rName);

                if (i == 0 || !(createdStructs[0 .. i].canFind(newMember)))
                {
                    structsBuf.put(GenDistStruct!newMember);
                    structsBuf.put("\n");
                    createdStructs[i] = newMember;
                    i++;
                }
            }
        }
    }
    return structsBuf.data;
}

mixin(GenDistStructs());

/// Calling distrib functions from Distribution structures
@system unittest
{
    import dstats.distrib;
    import std.math : approxEqual;

    assert(approxEqual(Normal.cdf(0.5, 0, 2), normalCDF(0.5, 0, 2)));
    assert(approxEqual(Beta.cdfr(0.5, 1, 1), betaCDFR(0.5, 1, 1)));
    assert(approxEqual(Cauchy.icdf(0.5, 0, 2), invCauchyCDF(0.5, 0, 2)));
    assert(approxEqual(Gamma.pdf(0.25, 1, 1), gammaPDF(0.25, 1, 1), ));
    assert(approxEqual(Hypergeometric.pmf(1, 2, 3, 5),
                       hypergeometricPMF(1, 2, 3, 5)));
    assert(approxEqual(Laplace.density(0.5, 0, 2), laplacePDF(0.5, 0, 2)));
    auto wei = Weibull.rand(1, 5);
}

@system unittest
{
    import dstats.distrib;
    import std.math : approxEqual;

    assert(approxEqual(betaCDF(0.5, 1, 1), Beta.cdf(0.5, 1, 1)));
    assert(approxEqual(betaCDFR(0.5, 1, 1), Beta.cdfr(0.5, 1, 1)));
    assert(approxEqual(invBetaCDF(0.5, 1, 1), Beta.icdf(0.5, 1, 1)));
    assert(approxEqual(betaPDF(0.5, 1, 1), Beta.pdf(0.5, 1, 1)));
    assert(approxEqual(betaPDF(0.5, 1, 1), Beta.density(0.5, 1, 1)));

    assert(approxEqual(binomialCDF(5, 10, 0.5), Binomial.cdf(5, 10, 0.5)));
    assert(approxEqual(binomialCDFR(5, 10, 0.5), Binomial.cdfr(5, 10, 0.5)));
    assert(approxEqual(invBinomialCDF(0.5, 10, 0.5),
                       Binomial.icdf(0.5, 10, 0.5)));
    assert(approxEqual(binomialPMF(5, 10, 0.5), Binomial.pmf(5, 10, 0.5)));
    assert(approxEqual(binomialPMF(5, 10, 0.5),
                       Binomial.density(5, 10, 0.5)));

    assert(approxEqual(cauchyCDF(0.5, 0, 2), Cauchy.cdf(0.5, 0, 2)));
    assert(approxEqual(cauchyCDFR(0.5, 0, 2), Cauchy.cdfr(0.5, 0, 2)));
    assert(approxEqual(invCauchyCDF(0.5, 0, 2), Cauchy.icdf(0.5, 0, 2)));
    assert(approxEqual(cauchyPDF(0.5, 0, 2), Cauchy.pdf(0.5, 0, 2)));
    assert(approxEqual(cauchyPDF(0.5, 0, 2), Cauchy.density(0.5, 0, 2)));

    assert(approxEqual(chiSquareCDF(0.5, 2), ChiSquare.cdf(0.5, 2)));
    assert(approxEqual(chiSquareCDFR(0.5, 2), ChiSquare.cdfr(0.5, 2)));
    assert(approxEqual(invChiSquareCDFR(2, 0.5), ChiSquare.icdfr(2, 0.5)));
    assert(approxEqual(chiSquarePDF(0.5, 2), ChiSquare.pdf(0.5, 2)));
    assert(approxEqual(chiSquarePDF(0.5, 2), ChiSquare.density(0.5, 2)));

    assert(approxEqual(dirichletPDF([0.25, 0.75], [0.25, 0.75]),
                       Dirichlet.pdf([0.25, 0.75], [0.25, 0.75])));

    assert(approxEqual(exponentialCDF(0.5, 1), Exponential.cdf(0.5, 1)));
    assert(approxEqual(exponentialCDFR(0.5, 1), Exponential.cdfr(0.5, 1)));
    assert(approxEqual(invExponentialCDF(0.5, 1), Exponential.icdf(0.5, 1)));
    assert(approxEqual(exponentialPDF(0.5, 1), Exponential.pdf(0.5, 1)));
    assert(approxEqual(exponentialPDF(0.5, 1), Exponential.density(0.5, 1)));

    assert(approxEqual(fisherCDF(0.25, 1, 1), Fisher.cdf(0.25, 1, 1)));
    assert(approxEqual(fisherCDFR(0.25, 1, 1), Fisher.cdfr(0.25, 1, 1)));
    assert(approxEqual(invFisherCDFR(1, 1, 0.25), Fisher.icdfr(1, 1, 0.25)));

    assert(approxEqual(gammaCDF(0.25, 1, 1), Gamma.cdf(0.25, 1, 1)));
    assert(approxEqual(gammaCDFR(0.25, 1, 1), Gamma.cdfr(0.25, 1, 1)));
    assert(approxEqual(invGammaCDF(0.25, 1, 1), Gamma.icdf(0.25, 1, 1)));
    assert(approxEqual(invGammaCDFR(0.25, 1, 1), Gamma.icdfr(0.25, 1, 1)));
    assert(approxEqual(gammaPDF(0.25, 1, 1), Gamma.pdf(0.25, 1, 1)));
    assert(approxEqual(gammaPDF(0.25, 1, 1), Gamma.density(0.25, 1, 1)));

    assert(approxEqual(hypergeometricCDF(1, 2, 3, 5),
                       Hypergeometric.cdf(1, 2, 3, 5)));
    assert(approxEqual(hypergeometricCDFR(1, 2, 3, 5),
                       Hypergeometric.cdfr(1, 2, 3, 5)));
    assert(approxEqual(hypergeometricPMF(1, 2, 3, 5),
                       Hypergeometric.pmf(1, 2, 3, 5)));
    assert(approxEqual(hypergeometricPMF(1, 2, 3, 5),
                       Hypergeometric.density(1, 2, 3, 5)));

    assert(approxEqual(laplaceCDF(0.5, 0, 2), Laplace.cdf(0.5, 0, 2)));
    assert(approxEqual(laplaceCDFR(0.5, 0, 2), Laplace.cdfr(0.5, 0, 2)));
    assert(approxEqual(invLaplaceCDF(0.5, 0, 2), Laplace.icdf(0.5, 0, 2)));
    assert(approxEqual(laplacePDF(0.5, 0, 2), Laplace.pdf(0.5, 0, 2)));
    assert(approxEqual(laplacePDF(0.5, 0, 2), Laplace.density(0.5, 0, 2)));

    assert(approxEqual(logisticCDF(0.5, 2, 5), Logistic.cdf(0.5, 2, 5)));

    assert(approxEqual(logNormalCDF(0.5, 0, 2), LogNormal.cdf(0.5, 0, 2)));
    assert(approxEqual(logNormalCDFR(0.5, 0, 2), LogNormal.cdfr(0.5, 0, 2)));
    assert(approxEqual(logNormalPDF(0.5, 0, 2), LogNormal.pdf(0.5, 0, 2)));
    assert(approxEqual(logNormalPDF(0.5, 0, 2),
                       LogNormal.density(0.5, 0, 2)));

    assert(approxEqual(negBinomCDF(5, 10, 0.5), NegBinom.cdf(5, 10, 0.5)));
    assert(approxEqual(negBinomCDFR(5, 10, 0.5), NegBinom.cdfr(5, 10, 0.5)));
    assert(approxEqual(invNegBinomCDF(0.5, 10, 0.5),
                       NegBinom.icdf(0.5, 10, 0.5)));
    assert(approxEqual(negBinomPMF(5, 10, 0.5), NegBinom.pmf(5, 10, 0.5)));
    assert(approxEqual(negBinomPMF(5, 10, 0.5),
                       NegBinom.density(5, 10, 0.5)));

    assert(approxEqual(normalCDF(0.5, 0, 2), Normal.cdf(0.5, 0, 2)));
    assert(approxEqual(normalCDFR(0.5, 0, 2), Normal.cdfr(0.5, 0, 2)));
    assert(approxEqual(invNormalCDF(0.5, 0, 2), Normal.icdf(0.5, 0, 2)));
    assert(approxEqual(normalPDF(0.5, 0, 2), Normal.pdf(0.5, 0, 2)));
    assert(approxEqual(normalPDF(0.5, 0, 2), Normal.density(0.5, 0, 2)));

    assert(approxEqual(poissonCDF(2, 4), Poisson.cdf(2, 4)));
    assert(approxEqual(poissonCDFR(2, 4), Poisson.cdfr(2, 4)));
    assert(approxEqual(invPoissonCDF(0.2, 4), Poisson.icdf(0.2, 4)));
    assert(approxEqual(poissonPMF(2, 4), Poisson.pmf(2, 4)));
    assert(approxEqual(poissonPMF(2, 4), Poisson.density(2, 4)));

    assert(approxEqual(rayleighCDF(0.5, 2), Rayleigh.cdf(0.5, 2)));

    assert(approxEqual(studentsTCDF(0.5, 10), StudentsT.cdf(0.5, 10)));
    assert(approxEqual(studentsTCDFR(0.5, 10), StudentsT.cdfr(0.5, 10)));
    assert(approxEqual(invStudentsTCDF(0.5, 10), StudentsT.icdf(0.5, 10)));
    assert(approxEqual(studentsTPDF(0.5, 10), StudentsT.pdf(0.5, 10)));
    assert(approxEqual(studentsTPDF(0.5, 10), StudentsT.density(0.5, 10)));

    assert(approxEqual(uniformCDF(0.5, 0, 2), Uniform.cdf(0.5, 0, 2)));
    assert(approxEqual(uniformCDFR(0.5, 0, 2), Uniform.cdfr(0.5, 0, 2)));
    assert(approxEqual(uniformPDF(0.5, 0, 2), Uniform.pdf(0.5, 0, 2)));
    assert(approxEqual(uniformPDF(0.5, 0, 2), Uniform.density(0.5, 0, 2)));

    assert(approxEqual(waldCDF(0.5, 1, 1), Wald.cdf(0.5, 1, 1)));

    assert(approxEqual(weibullCDF(0.5, 1, 5), Weibull.cdf(0.5, 1, 5)));
    assert(approxEqual(weibullCDFR(0.5, 1, 5), Weibull.cdfr(0.5, 1, 5)));
    assert(approxEqual(weibullPDF(0.5, 1, 5), Weibull.pdf(0.5, 1, 5)));
    assert(approxEqual(weibullPDF(0.5, 1, 5), Weibull.density(0.5, 1, 5)));
}

@system unittest
{
    import dstats.random;

    auto ber = Bernoulli.rand(0.5);
    auto bet = Beta.rand(1, 1);
    auto bin = Binomial.rand(10, 0.5);
    auto cau = Cauchy.rand(0, 2);
    auto chi = ChiSquare.rand(2);
    auto exp = Exponential.rand(1);
    auto fis = Fisher.rand(1, 1);
    auto gam = Gamma.rand(1, 1);
    auto geo = Geometric.rand(0.5);
    auto hyp = Hypergeometric.rand(2, 3, 4);
    auto lap = Laplace.rand(0, 2);
    auto log = Logistic.rand(2, 5);
    auto logn = LogNormal.rand(0, 2);
    auto neg = NegBinom.rand(10, 0.5);
    auto nor = Normal.rand(0, 2);
    auto poi = Poisson.rand(4);
    auto ray = Rayleigh.rand(2);
    auto stu = StudentsT.rand(10);
    auto wal = Wald.rand(1, 1);
    auto wei = Weibull.rand(1, 5);
}

private template GenDistFunc(string name)
{
    const char[] GenDistFunc =
        "///" ~ "\n" ~
        "auto " ~ name ~ "(alias T, U...)(U u)\n" ~
        "{\n" ~
        `    hasMemberCheck!(T, "` ~ name ~ `");` ~ "\n" ~
        "    return T." ~ name ~ "(u);\n" ~
        "}";
}

@system unittest
{
    auto val = GenDistFunc!"pdf";
    auto test =
        "///" ~ "\n" ~
        "auto " ~ "pdf" ~ "(alias T, U...)(U u)\n" ~
        "{\n" ~
        `    hasMemberCheck!(T, "` ~ "pdf" ~ `");` ~ "\n" ~
        "    return T." ~ "pdf" ~ "(u);\n" ~
        "}";
    assert(val == test);
}

mixin(GenDistFunc!("pdf"));
mixin(GenDistFunc!("pmf"));
mixin(GenDistFunc!("cdf"));
mixin(GenDistFunc!("cdfr"));
mixin(GenDistFunc!("icdf"));
mixin(GenDistFunc!("icdfr"));
mixin(GenDistFunc!("density"));
mixin(GenDistFunc!("rand"));

/// Using template functions to with behavior switching based on Distribution
@system unittest
{
    import dstats.distrib;
    import std.math : approxEqual;

    assert(approxEqual(cdf!Normal(0.5, 0, 2), normalCDF(0.5, 0, 2)));
    assert(approxEqual(cdfr!Beta(0.5, 1, 1), betaCDFR(0.5, 1, 1)));
    assert(approxEqual(icdf!Cauchy(0.5, 0, 2), invCauchyCDF(0.5, 0, 2)));
    assert(approxEqual(pdf!Gamma(0.25, 1, 1), gammaPDF(0.25, 1, 1)));
    assert(approxEqual(pmf!Hypergeometric(1, 2, 3, 5),
                       hypergeometricPMF(1, 2, 3, 5)));
    assert(approxEqual(density!Laplace(0.5, 0, 2), laplacePDF(0.5, 0, 2)));
    auto wei = rand!Weibull(1, 5);
}

@system unittest
{
    import std.math : approxEqual;

    assert(approxEqual(cdf!Beta(0.5, 1, 1), Beta.cdf(0.5, 1, 1)));
    assert(approxEqual(cdfr!Beta(0.5, 1, 1), Beta.cdfr(0.5, 1, 1)));
    assert(approxEqual(icdf!Beta(0.5, 1, 1), Beta.icdf(0.5, 1, 1)));
    assert(approxEqual(pdf!Beta(0.5, 1, 1), Beta.pdf(0.5, 1, 1)));
    assert(approxEqual(density!Beta(0.5, 1, 1), Beta.density(0.5, 1, 1)));

    assert(approxEqual(cdf!Binomial(5, 10, 0.5), Binomial.cdf(5, 10, 0.5)));
    assert(approxEqual(cdfr!Binomial(5, 10, 0.5), Binomial.cdfr(5, 10, 0.5)));
    assert(approxEqual(icdf!Binomial(0.5, 10, 0.5),
                       Binomial.icdf(0.5, 10, 0.5)));
    assert(approxEqual(pmf!Binomial(5, 10, 0.5), Binomial.pmf(5, 10, 0.5)));
    assert(approxEqual(density!Binomial(5, 10, 0.5),
                       Binomial.density(5, 10, 0.5)));

    assert(approxEqual(cdf!Cauchy(0.5, 0, 2), Cauchy.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr!Cauchy(0.5, 0, 2), Cauchy.cdfr(0.5, 0, 2)));
    assert(approxEqual(icdf!Cauchy(0.5, 0, 2), Cauchy.icdf(0.5, 0, 2)));
    assert(approxEqual(pdf!Cauchy(0.5, 0, 2), Cauchy.pdf(0.5, 0, 2)));
    assert(approxEqual(density!Cauchy(0.5, 0, 2), Cauchy.density(0.5, 0, 2)));

    assert(approxEqual(cdf!ChiSquare(0.5, 2), ChiSquare.cdf(0.5, 2)));
    assert(approxEqual(cdfr!ChiSquare(0.5, 2), ChiSquare.cdfr(0.5, 2)));
    assert(approxEqual(icdfr!ChiSquare(2, 0.5), ChiSquare.icdfr(2, 0.5)));
    assert(approxEqual(pdf!ChiSquare(0.5, 2), ChiSquare.pdf(0.5, 2)));
    assert(approxEqual(density!ChiSquare(0.5, 2), ChiSquare.density(0.5, 2)));

    assert(approxEqual(pdf!Dirichlet([0.25, 0.75], [0.25, 0.75]),
                       Dirichlet.pdf([0.25, 0.75], [0.25, 0.75])));

    assert(approxEqual(cdf!Exponential(0.5, 1), Exponential.cdf(0.5, 1)));
    assert(approxEqual(cdfr!Exponential(0.5, 1), Exponential.cdfr(0.5, 1)));
    assert(approxEqual(icdf!Exponential(0.5, 1), Exponential.icdf(0.5, 1)));
    assert(approxEqual(pdf!Exponential(0.5, 1), Exponential.pdf(0.5, 1)));
    assert(approxEqual(density!Exponential(0.5, 1),
                       Exponential.density(0.5, 1)));

    assert(approxEqual(cdf!Fisher(0.25, 1, 1), Fisher.cdf(0.25, 1, 1)));
    assert(approxEqual(cdfr!Fisher(0.25, 1, 1), Fisher.cdfr(0.25, 1, 1)));
    assert(approxEqual(icdfr!Fisher(1, 1, 0.25), Fisher.icdfr(1, 1, 0.25)));

    assert(approxEqual(cdf!Gamma(0.25, 1, 1), Gamma.cdf(0.25, 1, 1)));
    assert(approxEqual(cdfr!Gamma(0.25, 1, 1), Gamma.cdfr(0.25, 1, 1)));
    assert(approxEqual(icdf!Gamma(0.25, 1, 1), Gamma.icdf(0.25, 1, 1)));
    assert(approxEqual(icdfr!Gamma(0.25, 1, 1), Gamma.icdfr(0.25, 1, 1)));
    assert(approxEqual(pdf!Gamma(0.25, 1, 1), Gamma.pdf(0.25, 1, 1)));
    assert(approxEqual(density!Gamma(0.25, 1, 1), Gamma.density(0.25, 1, 1)));

    assert(approxEqual(cdf!Hypergeometric(1, 2, 3, 5),
                       Hypergeometric.cdf(1, 2, 3, 5)));
    assert(approxEqual(cdfr!Hypergeometric(1, 2, 3, 5),
                       Hypergeometric.cdfr(1, 2, 3, 5)));
    assert(approxEqual(pmf!Hypergeometric(1, 2, 3, 5),
                       Hypergeometric.pmf(1, 2, 3, 5)));
    assert(approxEqual(density!Hypergeometric(1, 2, 3, 5),
                       Hypergeometric.density(1, 2, 3, 5)));

    assert(approxEqual(cdf!Laplace(0.5, 0, 2), Laplace.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr!Laplace(0.5, 0, 2), Laplace.cdfr(0.5, 0, 2)));
    assert(approxEqual(icdf!Laplace(0.5, 0, 2), Laplace.icdf(0.5, 0, 2)));
    assert(approxEqual(pdf!Laplace(0.5, 0, 2), Laplace.pdf(0.5, 0, 2)));
    assert(approxEqual(density!Laplace(0.5, 0, 2), Laplace.density(0.5, 0, 2)));

    assert(approxEqual(cdf!Logistic(0.5, 2, 5), Logistic.cdf(0.5, 2, 5)));

    assert(approxEqual(cdf!LogNormal(0.5, 0, 2), LogNormal.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr!LogNormal(0.5, 0, 2), LogNormal.cdfr(0.5, 0, 2)));
    assert(approxEqual(pdf!LogNormal(0.5, 0, 2), LogNormal.pdf(0.5, 0, 2)));
    assert(approxEqual(density!LogNormal(0.5, 0, 2),
                       LogNormal.density(0.5, 0, 2)));

    assert(approxEqual(cdf!NegBinom(5, 10, 0.5), NegBinom.cdf(5, 10, 0.5)));
    assert(approxEqual(cdfr!NegBinom(5, 10, 0.5), NegBinom.cdfr(5, 10, 0.5)));
    assert(approxEqual(icdf!NegBinom(0.5, 10, 0.5),
                       NegBinom.icdf(0.5, 10, 0.5)));
    assert(approxEqual(pmf!NegBinom(5, 10, 0.5), NegBinom.pmf(5, 10, 0.5)));
    assert(approxEqual(density!NegBinom(5, 10, 0.5),
                       NegBinom.density(5, 10, 0.5)));

    assert(approxEqual(cdf!Normal(0.5, 0, 2), Normal.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr!Normal(0.5, 0, 2), Normal.cdfr(0.5, 0, 2)));
    assert(approxEqual(icdf!Normal(0.5, 0, 2), Normal.icdf(0.5, 0, 2)));
    assert(approxEqual(pdf!Normal(0.5, 0, 2), Normal.pdf(0.5, 0, 2)));
    assert(approxEqual(density!Normal(0.5, 0, 2), Normal.density(0.5, 0, 2)));

    assert(approxEqual(cdf!Poisson(2, 4), Poisson.cdf(2, 4)));
    assert(approxEqual(cdfr!Poisson(2, 4), Poisson.cdfr(2, 4)));
    assert(approxEqual(icdf!Poisson(0.2, 4), Poisson.icdf(0.2, 4)));
    assert(approxEqual(pmf!Poisson(2, 4), Poisson.pmf(2, 4)));
    assert(approxEqual(density!Poisson(2, 4), Poisson.density(2, 4)));

    assert(approxEqual(cdf!Rayleigh(0.5, 2), Rayleigh.cdf(0.5, 2)));

    assert(approxEqual(cdf!StudentsT(0.5, 10), StudentsT.cdf(0.5, 10)));
    assert(approxEqual(cdfr!StudentsT(0.5, 10), StudentsT.cdfr(0.5, 10)));
    assert(approxEqual(icdf!StudentsT(0.5, 10), StudentsT.icdf(0.5, 10)));
    assert(approxEqual(pdf!StudentsT(0.5, 10), StudentsT.pdf(0.5, 10)));
    assert(approxEqual(density!StudentsT(0.5, 10), StudentsT.density(0.5, 10)));

    assert(approxEqual(cdf!Uniform(0.5, 0, 2), Uniform.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr!Uniform(0.5, 0, 2), Uniform.cdfr(0.5, 0, 2)));
    assert(approxEqual(pdf!Uniform(0.5, 0, 2), Uniform.pdf(0.5, 0, 2)));
    assert(approxEqual(density!Uniform(0.5, 0, 2), Uniform.density(0.5, 0, 2)));

    assert(approxEqual(cdf!Wald(0.5, 1, 1), Wald.cdf(0.5, 1, 1)));

    assert(approxEqual(cdf!Weibull(0.5, 1, 5), Weibull.cdf(0.5, 1, 5)));
    assert(approxEqual(cdfr!Weibull(0.5, 1, 5), Weibull.cdfr(0.5, 1, 5)));
    assert(approxEqual(pdf!Weibull(0.5, 1, 5), Weibull.pdf(0.5, 1, 5)));
    assert(approxEqual(density!Weibull(0.5, 1, 5), Weibull.density(0.5, 1, 5)));
}

@system unittest
{
    auto ber = rand!Bernoulli(0.5);
    auto bet = rand!Beta(1, 1);
    auto bin = rand!Binomial(10, 0.5);
    auto cau = rand!Cauchy(0, 2);
    auto chi = rand!ChiSquare(2);
    auto exp = rand!Exponential(1);
    auto fis = rand!Fisher(1, 1);
    auto gam = rand!Gamma(1, 1);
    auto geo = rand!Geometric(0.5);
    auto hyp = rand!Hypergeometric(2, 3, 4);
    auto lap = rand!Laplace(0, 2);
    auto log = rand!Logistic(2, 5);
    auto logn = rand!LogNormal(0, 2);
    auto neg = rand!NegBinom(10, 0.5);
    auto nor = rand!Normal(0, 2);
    auto poi = rand!Poisson(4);
    auto ray = rand!Rayleigh(2);
    auto stu = rand!StudentsT(10);
    auto wal = rand!Wald(1, 1);
    auto wei = rand!Weibull(1, 5);
}

private template getParamType(alias T)
{
    import std.traits : isCallable, isFunction;

    static if (isFunction!(T) && isCallable!(T))
    {
        import std.traits : Parameters;
        alias getParamType = Parameters!(T)[1..$];
    }
    else
    {
        static assert ("SpecificDistribution only implements callable " ~
                       "functions, not templates.");
    }
}

unittest
{
    import dstats.distrib;
    alias ParametersType = getParamType!(betaPDF);
}

private template getParamTypeFull(alias T)
{
    static if (hasAnyDistrib!T)
    {
        static if (hasPDF!T)
            alias getParamTypeFull = getParamType!(T.pdf);
        else static if (hasPMF!T)
            alias getParamTypeFull = getParamType!(T.pmf);
        else static if (hasCDF!T)
            alias getParamTypeFull = getParamType!(T.cdf);
        else static if (hasCDFR!T)
            alias getParamTypeFull = getParamType!(T.cdfr);
        else static if (hasICDF!T)
            alias getParamTypeFull = getParamType!(T.icdf);
        else static if (hasICDFR!T)
            alias getParamTypeFull = getParamType!(T.icdfr);
        else
            static assert(0, "Should not be here");
    }
    else
    {
        static assert ("Type does not have any distribution functions.");
    }
}

private template getParamTypeRand(alias T)
{
    import std.traits : isCallable, isFunction;

    static if (isFunction!(T) && isCallable!(T))
    {
        import std.traits : Parameters;
        alias getParamTypeRand = Parameters!(T)[0..($ - 1)];
    }
    else
    {
        static assert ("SpecificDistribution only implements callable " ~
                       "functions, not templates.");
    }
}

private template aliasDist(Distribution, string s)
    if (isAggregateType!Distribution)
{
    static if (s == "cdf" && hasCDF!Distribution)
        alias aliasDist = Distribution.cdf;
    else static if (s == "pdf" && hasPDF!Distribution)
        alias aliasDist = Distribution.pdf;
    else static if (s == "pmf" && hasPMF!Distribution)
        alias aliasDist = Distribution.pmf;
    else static if (s == "density" && hasDensity!Distribution)
        alias aliasDist = Distribution.density;
    else static if (s == "icdf" && hasICDF!Distribution)
        alias aliasDist = Distribution.icdf;
    else static if (s == "cdfr" && hasCDFR!Distribution)
        alias aliasDist = Distribution.cdfr;
    else static if (s == "icdfr" && hasICDFR!Distribution)
        alias aliasDist = Distribution.icdfr;
    else
        static assert(0, "Member function " ~ s ~ " not implemented for " ~
                         "distribution: " ~ Distribution);
}

/**Struct wrapping Distribution behavior with parameters.
 *
 */
template SpecificDistribution(Distribution, RGen = Random)
    if (isAggregateType!Distribution && !hasAnyDistribTemplate!Distribution)
{
    struct SpecificDistribution
    {
        static if (hasAnyDistrib!Distribution)
        {
            alias ParametersType = getParamTypeFull!Distribution;
            ParametersType parameters;

            static if (hasRand!Distribution)
            {
                RGen gen;

                this(ParametersType x)
                {
                    import std.random : rndGen;

                    parameters = x;
                    gen = rndGen;
                }

                this(ParametersType x, RGen y)
                {
                    parameters = x;
                    gen = y;
                }
            }
            else
            {
                this(ParametersType x)
                {
                    parameters = x;
                }
            }
        }
        else static if (hasRand!Distribution)
        {
            alias ParametersType = getParamTypeRand!(Distribution.rand!RGen);

            ParametersType parameters;
            RGen gen;

            this(ParametersType x)
            {
                import std.random : rndGen;

                parameters = x;
                gen = rndGen;
            }

            this(ParametersType x, RGen y)
            {
                parameters = x;
                gen = y;
            }
        }
        else
        {
            static assert(0, Distribution ~ " is not a valid template " ~
                             "argument for " ~ __NAME__);
        }

        auto opDispatch(string s, U)(U x)
        {
            static if (s == "cdf" || s == "pdf" || s == "pmf" ||
                       s == "density" || s == "icdf" || s == "cdfr" ||
                       s == "icdfr")
            {
                alias f = aliasDist!(Distribution, s);
                alias fParameters = getParamType!f;
                static if (is(fParameters == ParametersType))
                {
                    return f(x, parameters);
                }
                else
                {
                    import std.conv : to;

                    fParameters fParam;
                    foreach(size_t i, ref e; fParam)
                        e = parameters[i].to!(fParameters[i]);
                    return f(x, fParam);
                }
            }
            else
            {
                static assert(0, "Member function " ~ s ~ " not implemented " ~
                              "for distribution: " ~ Distribution);
            }
        }

        auto opDispatch(string s)()
        {
            static if (s == "rand" && hasRand!Distribution)
            {
                alias rParameters = getParamTypeRand!(Distribution.rand!RGen);
                static if (is(rParameters == ParametersType))
                {
                    return Distribution.rand(parameters, gen);
                }
                else
                {
                    import std.conv : to;

                    rParameters rParam;
                    foreach(size_t i, ref e; rParam)
                        e = parameters[i].to!(rParameters[i]);
                    return Distribution.rand(rParam, gen);
                }
            }
            else
            {
                static assert(0, "Member function " ~ s ~ " not implemented " ~
                              "for distribution: " ~ Distriubtion);
            }
        }
    }
}

/**Handles the case where at least one of the distribution functions is a
 * template, a little bit of overlap with above.
 */
template SpecificDistribution(ParametersType, Distribution, RGen = Random)
    if (isAggregateType!Distribution && hasAnyDistribTemplate!Distribution)
{
    struct SpecificDistribution
    {
        static if (hasAnyDistrib!Distribution)
        {
            ParametersType parameters;

            static if (hasRand!Distribution)
            {
                RGen gen;

                this(auto ref ParametersType x) inout
                {
                    import std.random : rndGen;

                    parameters = x;
                    gen = rndGen;
                }

                this(ParametersType x, auto ref RGen y) inout
                {
                    parameters = x;
                    gen = y;
                }
            }
            else
            {
                this(ParametersType x)
                {
                    parameters = x;
                }
            }
        }
        else
        {
            static assert(0, Distribution ~ " is not a valid template " ~
                             "argument for " ~ __NAME__);
        }

        auto opDispatch(string s, U)(U x)
        {
            static if (s == "cdf" || s == "pdf" || s == "pmf" ||
                       s == "density" || s == "icdf" || s == "cdfr" ||
                       s == "icdfr")
            {
                alias f = aliasDist!(Distribution, s);
                static if (isTemplate!(f))
                {
                    return f!(U, ParametersType)(x, parameters);
                }
                else
                {
                    alias fParameters = getParamType!f;
                    static if (is(fParameters == ParametersType))
                    {
                        return f(x, parameters);
                    }
                    else
                    {
                        import std.conv : to;

                        fParameters fParam;
                        foreach(size_t i, ref e; fParam)
                            e = parameters[i].to!(fParameters[i]);
                        return f(x, fParam);
                    }
                }
            }
            else
            {
                static assert(0, "Member function " ~ s ~ " not implemented " ~
                              "for distribution: " ~ Distribution);
            }
        }

        auto opDispatch(string s)()
        {
            static if (s == "rand" && hasRand!Distribution)
            {
                static if (isTemplate!(Distribution.rand))
                {
                    return Distribution.rand!(ParametersType, RGen)
                                                              (parameters, gen);
                }
                else
                {
                    alias rParameters = getParamTypeRand!
                                                       (Distribution.rand!RGen);
                    static if (is(rParameters == ParametersType))
                    {
                        return Distribution.rand(parameters, gen);
                    }
                    else
                    {
                        import std.conv : to;

                        rParameters rParam;
                        foreach(size_t i, ref e; rParam)
                            e = parameters[i].to!(rParameters[i]);
                        return Distribution.rand(rParam, gen);
                    }
                }
            }
            else
            {
                static assert(0, "Member function " ~ s ~ " not implemented " ~
                              "for distribution: " ~ Distriubtion);
            }
        }
    }
}

/**Takes a distribution struct as a template argument, and
 * parameters as function arguments in the order that they appear in the
 * function declaration and returns a struct with those parameters as members
 * and that will forward any calls to the functions of the original distribution
 * to that distribution using the paramters supplied.
 *
 * Assumes the non-parameter argument is the first argument to the distribution
 * function (this is the case in nearly every function in dstats.distrb, but
 * invChiSquareCDFR and invFisherCDFR are currently exceptions).
 *
 * In addition, this function makes the assumption that the type signature of
 * the remaining parameters is the same across the different functions that
 * could be called. This causes a problem with invBinomialCDF currently as it
 * has a different type for the second parameter than other Binomial functions.
 *
 * Similar functionality to dstats.distrib.parametrize and
 * dstats.distrib.paramFunctor.
 */
template specificDistribution(Distribution, RGen = Random)
    if (isAggregateType!(Distribution))
{
    auto specificDistribution(Args...)(Args args)
    {
        static if (!hasAnyDistribTemplate!Distribution)
        {
            return SpecificDistribution!(Distribution, RGen)(args);
        }
        else
        {
            return SpecificDistribution!(Args, Distribution, RGen)(args);
        }
    }
}

/// Calling distrib functions from SpecificDistribution structures
@system unittest
{
    import std.math : approxEqual;

    auto stdNormal = specificDistribution!(Normal)(0, 1);

    assert(approxEqual(stdNormal.cdf(0.5), Normal.cdf(0.5, 0, 1)));
    assert(approxEqual(stdNormal.cdfr(0.5), Normal.cdfr(0.5, 0, 1)));
    assert(approxEqual(stdNormal.icdf(0.5), Normal.icdf(0.5, 0, 1)));
    assert(approxEqual(stdNormal.pdf(0.5), Normal.pdf(0.5, 0, 1)));
    assert(approxEqual(stdNormal.density(0.5), Normal.density(0.5, 0, 1)));
    auto rStdNorm = stdNormal.rand;
}

@system unittest
{
    import std.math : approxEqual;

    auto beta = specificDistribution!(Beta)(1, 1);
    assert(approxEqual(beta.cdf(0.5), Beta.cdf(0.5, 1, 1)));
    assert(approxEqual(beta.cdfr(0.5), Beta.cdfr(0.5, 1, 1)));
    assert(approxEqual(beta.icdf(0.5), Beta.icdf(0.5, 1, 1)));
    assert(approxEqual(beta.pdf(0.5), Beta.pdf(0.5, 1, 1)));
    assert(approxEqual(beta.density(0.5), Beta.density(0.5, 1, 1)));

    auto bino = specificDistribution!(Binomial)(10, 0.5);
    assert(approxEqual(bino.cdf(5), Binomial.cdf(5, 10, 0.5)));
    assert(approxEqual(bino.cdfr(5), Binomial.cdfr(5, 10, 0.5)));
    assert(approxEqual(bino.icdf(0.5), Binomial.icdf(0.5, 10, 0.5)));
    assert(approxEqual(bino.pmf(5), Binomial.pmf(5, 10, 0.5)));
    assert(approxEqual(bino.density(5), Binomial.density(5, 10, 0.5)));

    auto cauc = specificDistribution!(Cauchy)(0, 2);
    assert(approxEqual(cauc.cdf(0.5), Cauchy.cdf(0.5, 0, 2)));
    assert(approxEqual(cauc.cdfr(0.5), Cauchy.cdfr(0.5, 0, 2)));
    assert(approxEqual(cauc.icdf(0.5), Cauchy.icdf(0.5, 0, 2)));
    assert(approxEqual(cauc.pdf(0.5), Cauchy.pdf(0.5, 0, 2)));
    assert(approxEqual(cauc.density(0.5), Cauchy.density(0.5, 0, 2)));

    auto chis = specificDistribution!(ChiSquare)(2);
    assert(approxEqual(chis.cdf(0.5), ChiSquare.cdf(0.5, 2)));
    assert(approxEqual(chis.cdfr(0.5), ChiSquare.cdfr(0.5, 2)));
    //assert(approxEqual(chis.icdfr(0.5), ChiSquare.icdfr(2, 0.5)));
    assert(approxEqual(chis.pdf(0.5), ChiSquare.pdf(0.5, 2)));
    assert(approxEqual(chis.density(0.5), ChiSquare.density(0.5, 2)));

    auto diri = specificDistribution!(Dirichlet)([0.25, 0.75]);
    assert(approxEqual(diri.pdf([0.25, 0.75]),
                       Dirichlet.pdf([0.25, 0.75], [0.25, 0.75])));

    auto expo = specificDistribution!(Exponential)(1);
    assert(approxEqual(expo.cdf(0.5), Exponential.cdf(0.5, 1)));
    assert(approxEqual(expo.cdfr(0.5), Exponential.cdfr(0.5, 1)));
    assert(approxEqual(expo.icdf(0.5), Exponential.icdf(0.5, 1)));
    assert(approxEqual(expo.pdf(0.5), Exponential.pdf(0.5, 1)));
    assert(approxEqual(expo.density(0.5), Exponential.density(0.5, 1)));

    auto fish = specificDistribution!(Fisher)(1, 1);
    assert(approxEqual(fish.cdf(0.25), Fisher.cdf(0.25, 1, 1)));
    assert(approxEqual(fish.cdfr(0.25), Fisher.cdfr(0.25, 1, 1)));
    //assert(approxEqual(fish.icdfr(0.5), Fisher.icdfr(1, 1, 0.5)));

    auto gamm = specificDistribution!(Gamma)(1, 1);
    assert(approxEqual(gamm.cdf(0.25), Gamma.cdf(0.25, 1, 1)));
    assert(approxEqual(gamm.cdfr(0.25), Gamma.cdfr(0.25, 1, 1)));
    assert(approxEqual(gamm.icdf(0.25), Gamma.icdf(0.25, 1, 1)));
    assert(approxEqual(gamm.icdfr(0.25), Gamma.icdfr(0.25, 1, 1)));
    assert(approxEqual(gamm.pdf(0.25), Gamma.pdf(0.25, 1, 1)));
    assert(approxEqual(gamm.density(0.25), Gamma.density(0.25, 1, 1)));

    auto hype = specificDistribution!(Hypergeometric)(2, 3, 4);
    assert(approxEqual(hype.cdf(1), Hypergeometric.cdf(1, 2, 3, 4)));
    assert(approxEqual(hype.cdfr(1), Hypergeometric.cdfr(1, 2, 3, 4)));
    assert(approxEqual(hype.pmf(1), Hypergeometric.pmf(1, 2, 3, 4)));
    assert(approxEqual(hype.density(1), Hypergeometric.density(1, 2, 3, 4)));

    auto lapl = specificDistribution!(Laplace)(0, 2);
    assert(approxEqual(lapl.cdf(0.5), Laplace.cdf(0.5, 0, 2)));
    assert(approxEqual(lapl.cdfr(0.5), Laplace.cdfr(0.5, 0, 2)));
    assert(approxEqual(lapl.icdf(0.5), Laplace.icdf(0.5, 0, 2)));
    assert(approxEqual(lapl.pdf(0.5), Laplace.pdf(0.5, 0, 2)));
    assert(approxEqual(lapl.density(0.5), Laplace.density(0.5, 0, 2)));

    auto logi = specificDistribution!(Logistic)(2, 5);
    assert(approxEqual(logi.cdf(0.5), Logistic.cdf(0.5, 2, 5)));

    auto logn = specificDistribution!(LogNormal)(0, 2);
    assert(approxEqual(logn.cdf(0.5), LogNormal.cdf(0.5, 0, 2)));
    assert(approxEqual(logn.cdfr(0.5), LogNormal.cdfr(0.5, 0, 2)));
    assert(approxEqual(logn.pdf(0.5), LogNormal.pdf(0.5, 0, 2)));
    assert(approxEqual(logn.density(0.5), LogNormal.density(0.5, 0, 2)));

    auto negb = specificDistribution!(NegBinom)(10, 0.5);
    assert(approxEqual(negb.cdf(5), NegBinom.cdf(5, 10, 0.5)));
    assert(approxEqual(negb.cdfr(5), NegBinom.cdfr(5, 10, 0.5)));
    assert(approxEqual(negb.icdf(0.5), NegBinom.icdf(0.5, 10, 0.5)));
    assert(approxEqual(negb.pmf(5), NegBinom.pmf(5, 10, 0.5)));
    assert(approxEqual(negb.density(5), NegBinom.density(5, 10, 0.5)));

    auto norm = specificDistribution!(Normal)(0, 2);
    assert(approxEqual(norm.cdf(0.5), Normal.cdf(0.5, 0, 2)));
    assert(approxEqual(norm.cdfr(0.5), Normal.cdfr(0.5, 0, 2)));
    assert(approxEqual(norm.icdf(0.5), Normal.icdf(0.5, 0, 2)));
    assert(approxEqual(norm.pdf(0.5), Normal.pdf(0.5, 0, 2)));
    assert(approxEqual(norm.density(0.5), Normal.density(0.5, 0, 2)));

    auto pois = specificDistribution!(Poisson)(4);
    assert(approxEqual(pois.cdf(2), Poisson.cdf(2, 4)));
    assert(approxEqual(pois.cdfr(2), Poisson.cdfr(2, 4)));
    assert(approxEqual(pois.icdf(0.2), Poisson.icdf(0.2, 4)));
    assert(approxEqual(pois.pmf(2), Poisson.pmf(2, 4)));
    assert(approxEqual(pois.density(2), Poisson.density(2, 4)));

    auto rayl = specificDistribution!(Rayleigh)(2);
    assert(approxEqual(rayl.cdf(0.5), Rayleigh.cdf(0.5, 2)));

    auto stud = specificDistribution!(StudentsT)(10);
    assert(approxEqual(stud.cdf(0.5), StudentsT.cdf(0.5, 10)));
    assert(approxEqual(stud.cdfr(0.5), StudentsT.cdfr(0.5, 10)));
    assert(approxEqual(stud.icdf(0.5), StudentsT.icdf(0.5, 10)));
    assert(approxEqual(stud.pdf(0.5), StudentsT.pdf(0.5, 10)));
    assert(approxEqual(stud.density(0.5), StudentsT.density(0.5, 10)));

    auto unif = specificDistribution!(Uniform)(0, 2);
    assert(approxEqual(unif.cdf(0.5), Uniform.cdf(0.5, 0, 2)));
    assert(approxEqual(unif.cdfr(0.5), Uniform.cdfr(0.5, 0, 2)));
    assert(approxEqual(unif.pdf(0.5), Uniform.pdf(0.5, 0, 2)));
    assert(approxEqual(unif.density(0.5), Uniform.density(0.5, 0, 2)));

    auto wal = specificDistribution!(Wald)(1, 1);
    assert(approxEqual(wal.cdf(0.5), Wald.cdf(0.5, 1, 1)));

    auto weib = specificDistribution!(Weibull)(1, 5);
    assert(approxEqual(weib.cdf(0.5), Weibull.cdf(0.5, 1, 5)));
    assert(approxEqual(weib.cdfr(0.5), Weibull.cdfr(0.5, 1, 5)));
    assert(approxEqual(weib.pdf(0.5), Weibull.pdf(0.5, 1, 5)));
    assert(approxEqual(weib.density(0.5), Weibull.density(0.5, 1, 5)));
}

@system unittest
{
    auto bern = specificDistribution!(Bernoulli)(0.5);
    auto ber = Bernoulli.rand;

    auto beta = specificDistribution!(Beta)(1, 1);
    auto rbeta = beta.rand;

    auto bino = specificDistribution!(Binomial)(10, 0.5);
    auto rbino = bino.rand;

    auto cauc = specificDistribution!(Cauchy)(0, 2);
    auto rcauc = cauc.rand;

    auto chis = specificDistribution!(ChiSquare)(2);
    auto rchis = chis.rand;

    auto expo = specificDistribution!(Exponential)(1);
    auto rexpo = expo.rand;

    auto fish = specificDistribution!(Fisher)(1, 1);
    auto rfish = fish.rand;

    auto gamm = specificDistribution!(Gamma)(1, 1);
    auto rgamm = gamm.rand;

    auto hype = specificDistribution!(Hypergeometric)(2, 3, 4);
    auto rhype = hype.rand;

    auto lapl = specificDistribution!(Laplace)(0, 2);
    auto rlapl = lapl.rand;

    auto logi = specificDistribution!(Logistic)(2, 5);
    auto rlogi = logi.rand;

    auto logn = specificDistribution!(LogNormal)(0, 2);
    auto rlogn = logn.rand;

    auto negb = specificDistribution!(NegBinom)(10, 0.5);
    auto rnegb = negb.rand;

    auto norm = specificDistribution!(Normal)(0, 2);
    auto rnorm = norm.rand;

    auto pois = specificDistribution!(Poisson)(4);
    auto rpois = pois.rand;

    auto rayl = specificDistribution!(Rayleigh)(2);
    auto rrayl = rayl.rand;

    auto stud = specificDistribution!(StudentsT)(10);
    auto rstud = stud.rand;

    auto wal = specificDistribution!(Wald)(1, 1);
    auto rwal = wal.rand;

    auto weib = specificDistribution!(Weibull)(1, 5);
    auto rweib = weib.rand;
}

@system unittest
{
    import std.random : Random, rndGen;

    Random gen = rndGen;
    auto norm = specificDistribution!(Normal)(0, 2, gen);
    auto rnorm = norm.rand;

    auto bern = specificDistribution!(Bernoulli)(0.5, gen);
    auto ber = Bernoulli.rand(0.5);
}

private template GenDistFuncSpecificDistribution(string name)
{
    const char[] GenDistFuncSpecificDistribution =
        "///" ~ "\n" ~
        "auto " ~ name ~ "(T, U)(T x, U y)\n" ~
        "    if (isSpecificDistribution!U)\n" ~
        "{\n" ~
        "    return y." ~ name ~ "(x);\n" ~
        "}";
}

@system unittest
{
    auto val = GenDistFuncSpecificDistribution!"pdf";
    auto test =
        "///" ~ "\n" ~
        "auto " ~ "pdf" ~ "(T, U)(T x, U y)\n" ~
        "    if (isSpecificDistribution!U)\n" ~
        "{\n" ~
        "    return y.pdf(x);\n" ~
        "}";
    assert(val == test);
}

mixin(GenDistFuncSpecificDistribution!("cdf"));
mixin(GenDistFuncSpecificDistribution!("pdf"));
mixin(GenDistFuncSpecificDistribution!("pmf"));
mixin(GenDistFuncSpecificDistribution!("cdfr"));
mixin(GenDistFuncSpecificDistribution!("icdf"));
mixin(GenDistFuncSpecificDistribution!("icdfr"));
mixin(GenDistFuncSpecificDistribution!("density"));

auto rand(T)(T x)
    if (isSpecificDistribution!T)
{
    return x.rand;
}

/// Distribution functions that can call with SpecificDistribution parameters
@system unittest
{
    import std.math : approxEqual;

    auto stdNormal = specificDistribution!(Normal)(0, 1);

    assert(approxEqual(cdf(0.5, stdNormal), Normal.cdf(0.5, 0, 1)));
    assert(approxEqual(cdfr(0.5, stdNormal), Normal.cdfr(0.5, 0, 1)));
    assert(approxEqual(icdf(0.5, stdNormal), Normal.icdf(0.5, 0, 1)));
    assert(approxEqual(pdf(0.5, stdNormal), Normal.pdf(0.5, 0, 1)));
    assert(approxEqual(density(0.5, stdNormal), Normal.density(0.5, 0, 1)));
    auto rStdNorm = rand(stdNormal);
}

@system unittest
{
    import std.math : approxEqual;

    auto beta = specificDistribution!(Beta)(1, 1);
    assert(approxEqual(cdf(0.5, beta), Beta.cdf(0.5, 1, 1)));
    assert(approxEqual(cdfr(0.5, beta), Beta.cdfr(0.5, 1, 1)));
    assert(approxEqual(icdf(0.5, beta), Beta.icdf(0.5, 1, 1)));
    assert(approxEqual(pdf(0.5, beta), Beta.pdf(0.5, 1, 1)));
    assert(approxEqual(density(0.5, beta), Beta.density(0.5, 1, 1)));

    auto bino = specificDistribution!(Binomial)(10, 0.5);
    assert(approxEqual(cdf(5, bino), Binomial.cdf(5, 10, 0.5)));
    assert(approxEqual(cdfr(5, bino), Binomial.cdfr(5, 10, 0.5)));
    assert(approxEqual(icdf(0.5, bino), Binomial.icdf(0.5, 10, 0.5)));
    assert(approxEqual(pmf(5, bino), Binomial.pmf(5, 10, 0.5)));
    assert(approxEqual(density(5, bino), Binomial.density(5, 10, 0.5)));

    auto cauc = specificDistribution!(Cauchy)(0, 2);
    assert(approxEqual(cdf(0.5, cauc), Cauchy.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr(0.5, cauc), Cauchy.cdfr(0.5, 0, 2)));
    assert(approxEqual(icdf(0.5, cauc), Cauchy.icdf(0.5, 0, 2)));
    assert(approxEqual(pdf(0.5, cauc), Cauchy.pdf(0.5, 0, 2)));
    assert(approxEqual(density(0.5, cauc), Cauchy.density(0.5, 0, 2)));

    auto chis = specificDistribution!(ChiSquare)(2);
    assert(approxEqual(cdf(0.5, chis), ChiSquare.cdf(0.5, 2)));
    assert(approxEqual(cdfr(0.5, chis), ChiSquare.cdfr(0.5, 2)));
    //assert(approxEqual(icdfr(0.5, chis), ChiSquare.icdfr(2, 0.5)));
    assert(approxEqual(pdf(0.5, chis), ChiSquare.pdf(0.5, 2)));
    assert(approxEqual(density(0.5, chis), ChiSquare.density(0.5, 2)));

    auto diri = specificDistribution!(Dirichlet)([0.25, 0.75]);
    assert(approxEqual(pdf([0.25, 0.75], diri),
                       Dirichlet.pdf([0.25, 0.75], [0.25, 0.75])));

    auto expo = specificDistribution!(Exponential)(1);
    assert(approxEqual(cdf(0.5, expo), Exponential.cdf(0.5, 1)));
    assert(approxEqual(cdfr(0.5, expo), Exponential.cdfr(0.5, 1)));
    assert(approxEqual(icdf(0.5, expo), Exponential.icdf(0.5, 1)));
    assert(approxEqual(pdf(0.5, expo), Exponential.pdf(0.5, 1)));
    assert(approxEqual(density(0.5, expo), Exponential.density(0.5, 1)));

    auto fish = specificDistribution!(Fisher)(1, 1);
    assert(approxEqual(cdf(0.25, fish), Fisher.cdf(0.25, 1, 1)));
    assert(approxEqual(cdfr(0.25, fish), Fisher.cdfr(0.25, 1, 1)));
    //assert(approxEqual(icdfr(0.5, fish), Fisher.icdfr(1, 1, 0.5)));

    auto gamm = specificDistribution!(Gamma)(1, 1);
    assert(approxEqual(cdf(0.25, gamm), Gamma.cdf(0.25, 1, 1)));
    assert(approxEqual(cdfr(0.25, gamm), Gamma.cdfr(0.25, 1, 1)));
    assert(approxEqual(icdf(0.25, gamm), Gamma.icdf(0.25, 1, 1)));
    assert(approxEqual(icdfr(0.25, gamm), Gamma.icdfr(0.25, 1, 1)));
    assert(approxEqual(pdf(0.25, gamm), Gamma.pdf(0.25, 1, 1)));
    assert(approxEqual(density(0.25, gamm), Gamma.density(0.25, 1, 1)));

    auto hype = specificDistribution!(Hypergeometric)(2, 3, 4);
    assert(approxEqual(cdf(1, hype), Hypergeometric.cdf(1, 2, 3, 4)));
    assert(approxEqual(cdfr(1, hype), Hypergeometric.cdfr(1, 2, 3, 4)));
    assert(approxEqual(pmf(1, hype), Hypergeometric.pmf(1, 2, 3, 4)));
    assert(approxEqual(density(1, hype), Hypergeometric.density(1, 2, 3, 4)));

    auto lapl = specificDistribution!(Laplace)(0, 2);
    assert(approxEqual(cdf(0.5, lapl), Laplace.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr(0.5, lapl), Laplace.cdfr(0.5, 0, 2)));
    assert(approxEqual(icdf(0.5, lapl), Laplace.icdf(0.5, 0, 2)));
    assert(approxEqual(pdf(0.5, lapl), Laplace.pdf(0.5, 0, 2)));
    assert(approxEqual(density(0.5, lapl), Laplace.density(0.5, 0, 2)));

    auto logi = specificDistribution!(Logistic)(2, 5);
    assert(approxEqual(cdf(0.5, logi), Logistic.cdf(0.5, 2, 5)));

    auto logn = specificDistribution!(LogNormal)(0, 2);
    assert(approxEqual(cdf(0.5, logn), LogNormal.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr(0.5, logn), LogNormal.cdfr(0.5, 0, 2)));
    assert(approxEqual(pdf(0.5, logn), LogNormal.pdf(0.5, 0, 2)));
    assert(approxEqual(density(0.5, logn), LogNormal.density(0.5, 0, 2)));

    auto negb = specificDistribution!(NegBinom)(10, 0.5);
    assert(approxEqual(cdf(5, negb), NegBinom.cdf(5, 10, 0.5)));
    assert(approxEqual(cdfr(5, negb), NegBinom.cdfr(5, 10, 0.5)));
    assert(approxEqual(icdf(0.5, negb), NegBinom.icdf(0.5, 10, 0.5)));
    assert(approxEqual(pmf(5, negb), NegBinom.pmf(5, 10, 0.5)));
    assert(approxEqual(density(5, negb), NegBinom.density(5, 10, 0.5)));

    auto norm = specificDistribution!(Normal)(0, 2);
    assert(approxEqual(cdf(0.5, norm), Normal.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr(0.5, norm), Normal.cdfr(0.5, 0, 2)));
    assert(approxEqual(icdf(0.5, norm), Normal.icdf(0.5, 0, 2)));
    assert(approxEqual(pdf(0.5, norm), Normal.pdf(0.5, 0, 2)));
    assert(approxEqual(density(0.5, norm), Normal.density(0.5, 0, 2)));

    auto pois = specificDistribution!(Poisson)(4);
    assert(approxEqual(cdf(2, pois), Poisson.cdf(2, 4)));
    assert(approxEqual(cdfr(2, pois), Poisson.cdfr(2, 4)));
    assert(approxEqual(icdf(0.2, pois), Poisson.icdf(0.2, 4)));
    assert(approxEqual(pmf(2, pois), Poisson.pmf(2, 4)));
    assert(approxEqual(density(2, pois), Poisson.density(2, 4)));

    auto rayl = specificDistribution!(Rayleigh)(2);
    assert(approxEqual(cdf(0.5, rayl), Rayleigh.cdf(0.5, 2)));

    auto stud = specificDistribution!(StudentsT)(10);
    assert(approxEqual(cdf(0.5, stud), StudentsT.cdf(0.5, 10)));
    assert(approxEqual(cdfr(0.5, stud), StudentsT.cdfr(0.5, 10)));
    assert(approxEqual(icdf(0.5, stud), StudentsT.icdf(0.5, 10)));
    assert(approxEqual(pdf(0.5, stud), StudentsT.pdf(0.5, 10)));
    assert(approxEqual(density(0.5, stud), StudentsT.density(0.5, 10)));

    auto unif = specificDistribution!(Uniform)(0, 2);
    assert(approxEqual(cdf(0.5, unif), Uniform.cdf(0.5, 0, 2)));
    assert(approxEqual(cdfr(0.5, unif), Uniform.cdfr(0.5, 0, 2)));
    assert(approxEqual(pdf(0.5, unif), Uniform.pdf(0.5, 0, 2)));
    assert(approxEqual(density(0.5, unif), Uniform.density(0.5, 0, 2)));

    auto wal = specificDistribution!(Wald)(1, 1);
    assert(approxEqual(cdf(0.5, wal), Wald.cdf(0.5, 1, 1)));

    auto weib = specificDistribution!(Weibull)(1, 5);
    assert(approxEqual(cdf(0.5, weib), Weibull.cdf(0.5, 1, 5)));
    assert(approxEqual(cdfr(0.5, weib), Weibull.cdfr(0.5, 1, 5)));
    assert(approxEqual(pdf(0.5, weib), Weibull.pdf(0.5, 1, 5)));
    assert(approxEqual(density(0.5, weib), Weibull.density(0.5, 1, 5)));
}

@system unittest
{
    auto bern = specificDistribution!(Bernoulli)(0.5);
    auto ber = rand(bern);

    auto beta = specificDistribution!(Beta)(1, 1);
    auto rbeta = rand(beta);

    auto bino = specificDistribution!(Binomial)(10, 0.5);
    auto rbino = rand(bino);

    auto cauc = specificDistribution!(Cauchy)(0, 2);
    auto rcauc = rand(cauc);

    auto chis = specificDistribution!(ChiSquare)(2);
    auto rchis = rand(chis);

    auto expo = specificDistribution!(Exponential)(1);
    auto rexpo = rand(expo);

    auto fish = specificDistribution!(Fisher)(1, 1);
    auto rfish = rand(fish);

    auto gamm = specificDistribution!(Gamma)(1, 1);
    auto rgamm = rand(gamm);

    auto hype = specificDistribution!(Hypergeometric)(2, 3, 4);
    auto rhype = rand(hype);

    auto lapl = specificDistribution!(Laplace)(0, 2);
    auto rlapl = rand(lapl);

    auto logi = specificDistribution!(Logistic)(2, 5);
    auto rlogi = rand(logi);

    auto logn = specificDistribution!(LogNormal)(0, 2);
    auto rlogn = rand(logn);

    auto negb = specificDistribution!(NegBinom)(10, 0.5);
    auto rnegb = rand(negb);

    auto norm = specificDistribution!(Normal)(0, 2);
    auto rnorm = rand(norm);

    auto pois = specificDistribution!(Poisson)(4);
    auto rpois = rand(pois);

    auto rayl = specificDistribution!(Rayleigh)(2);
    auto rrayl = rand(rayl);

    auto stud = specificDistribution!(StudentsT)(10);
    auto rstud = rand(stud);

    auto wal = specificDistribution!(Wald)(1, 1);
    auto rwal = rand(wal);

    auto weib = specificDistribution!(Weibull)(1, 5);
    auto rweib = rand(weib);
}

private template isSpecificDistribution(T)
{
    static if (is(T : SpecificDistribution!(Distribution, RGen),
                      Distribution, RGen))
    {
        enum bool isSpecificDistribution = true;
    }
    else static if (is(T : SpecificDistribution!(ParametersType, Distribution,
                       RGen), ParametersType, Distribution, RGen))
    {
        enum bool isSpecificDistribution = true;
    }
    else
    {
        enum bool isSpecificDistribution = false;
    }
}

@safe unittest
{
    auto weib = specificDistribution!(Weibull)(1, 5);

    static assert(isSpecificDistribution!(typeof(weib)));
}

private template getDistributionOf(T)
    if (isSpecificDistribution!T)
{
    static if (is(T : SpecificDistribution!(Distribution, RGen),
                      Distribution, RGen))
    {
        alias getDistributionOf = Distribution;
    }
    else static if (is(T : SpecificDistribution!(ParametersType, Distribution,
                       RGen), ParametersType, Distribution, RGen))
    {
        alias getDistributionOf = Distribution;
    }
    else
    {
        static assert(0, "should not be here");
    }
}

@safe unittest
{
    auto weib = specificDistribution!(Weibull)(1, 5);

    static assert(is(getDistributionOf!(typeof(weib)) == Weibull));
}


/**Convenience function to allow one-statement creation of arrays of random
 * numbers.
 *
 * For more information refer to dstats.random.randArray.
 */
template randArray(Distribution)
    if (isAggregateType!(Distribution) &&
        !isSpecificDistribution!Distribution)
{
    auto randArray(Args...)(size_t N, auto ref Args args)
    {
        static import dstats.random;
        alias typeof(Distribution.rand(args)) R;
        return dstats.random.randArray!(R, Distribution.rand, Args)(N, args);
    }
}

///
template randArray(R, Distribution)
    if (isAggregateType!(Distribution) &&
        !isSpecificDistribution!Distribution)
{
    R[] randArray(Args...)(size_t N, auto ref Args args)
    {
        static import dstats.random;
        return dstats.random.randArray!(R, Distribution.rand, Args)(N, args);
    }
}

///
auto randArray(T : SpecificDistribution!(Distribution, RGen),
                   Distribution, RGen)(size_t N, T x)
    if (isSpecificDistribution!T)
{
    static import dstats.random;
    alias typeof(Distribution.rand(x.parameters, x.gen)) R;
    return dstats.random.randArray!(R, Distribution.rand)
                                                       (N, x.parameters, x.gen);
}

///
R[] randArray(R, T : SpecificDistribution!(Distribution, RGen),
                   Distribution, RGen)(size_t N, T x)
    if (isSpecificDistribution!T)
{
    static import dstats.random;
    return dstats.random.randArray!(R, Distribution.rand)
                                                       (N, x.parameters, x.gen);
}

///
auto randArray(T : SpecificDistribution!(ParametersType, Distribution, RGen),
                   ParametersType, Distribution, RGen)(size_t N, T x)
    if (isSpecificDistribution!T)
{
    static import dstats.random;
    alias typeof(Distribution.rand(x.parameters, x.gen)) R;
    return dstats.random.randArray!(R, Distribution.rand)
                                                       (N, x.parameters, x.gen);
}

///
R[] randArray(R, T : SpecificDistribution!(ParametersType, Distribution, RGen),
                     ParametersType, Distribution, RGen)(size_t N, T x)
    if (isSpecificDistribution!T)
{
    static import dstats.random;
    return dstats.random.randArray!(R, Distribution.rand)
                                                       (N, x.parameters, x.gen);
}

///
@system unittest
{
    enum size_t N = 5;

    auto x = randArray!Normal(N, 0, 1);
}

///
@system unittest
{
    enum size_t N = 5;

    auto stdNormal = specificDistribution!(Normal)(0, 1);
    auto x = randArray(N, stdNormal);
    float[] y = randArray!(float)(N, stdNormal);
}

@system unittest
{
    enum size_t N = 5;

    auto ber = randArray!Bernoulli(N, 0.5);
    auto bet = randArray!Beta(N, 1, 1);
    auto bin = randArray!Binomial(N, 10, 0.5);
    auto cau = randArray!Cauchy(N, 0, 2);
    auto chi = randArray!ChiSquare(N, 2);
    auto exp = randArray!Exponential(N, 1);
    auto fis = randArray!Fisher(N, 1, 1);
    auto gam = randArray!Gamma(N, 1, 1);
    auto geo = randArray!Geometric(N, 0.5);
    auto hyp = randArray!Hypergeometric(N, 2, 3, 5);
    auto lap = randArray!Laplace(N, 0, 2);
    auto log = randArray!Logistic(N, 2, 5);
    auto logn = randArray!LogNormal(N, 0, 2);
    auto neg = randArray!NegBinom(N, 10, 0.5);
    auto nor = randArray!Normal(N, 0, 2);
    auto poi = randArray!Poisson(N, 4);
    auto ray = randArray!Rayleigh(N, 2);
    auto stu = randArray!StudentsT(N, 10);
    auto wal = randArray!Wald(N, 1, 1);
    auto wei = randArray!Weibull(N, 1, 5);
}

@system unittest
{
    enum size_t N = 5;

    auto ber = randArray!(float, Bernoulli)(N, 0.5);
    auto bet = randArray!(float, Beta)(N, 1, 1);
    auto bin = randArray!(float, Binomial)(N, 10, 0.5);
    auto cau = randArray!(float, Cauchy)(N, 0, 2);
    auto chi = randArray!(float, ChiSquare)(N, 2);
    auto exp = randArray!(float, Exponential)(N, 1);
    auto fis = randArray!(float, Fisher)(N, 1, 1);
    auto gam = randArray!(float, Gamma)(N, 1, 1);
    auto geo = randArray!(float, Geometric)(N, 0.5);
    auto hyp = randArray!(float, Hypergeometric)(N, 2, 3, 5);
    auto lap = randArray!(float, Laplace)(N, 0, 2);
    auto log = randArray!(float, Logistic)(N, 2, 5);
    auto logn = randArray!(float, LogNormal)(N, 0, 2);
    auto neg = randArray!(float, NegBinom)(N, 10, 0.5);
    auto nor = randArray!(float, Normal)(N, 0, 2);
    auto poi = randArray!(float, Poisson)(N, 4);
    auto ray = randArray!(float, Rayleigh)(N, 2);
    auto stu = randArray!(float, StudentsT)(N, 10);
    auto wal = randArray!(float, Wald)(N, 1, 1);
    auto wei = randArray!(float, Weibull)(N, 1, 5);
}

/**Turn a random number generator function into an infinite range.
 * Params is a tuple of the distribution parameters.  This is specified
 * in the same order as when calling the function directly.
 *
 * For more information refer to dstats.random.randRange.
 *
 */
template randRange(Distribution)
    if (isAggregateType!Distribution &&
        !isSpecificDistribution!Distribution)
{
    auto randRange(T...)(T params)
    {
        static import dstats.random;
        return dstats.random.randRange!(Distribution.rand)(params);
    }
}

///
auto randRange(T : SpecificDistribution!(Distribution, RGen),
                   Distribution, RGen)(T x)
    if (isSpecificDistribution!T)
{
    static import dstats.random;
    return dstats.random.randRange!(Distribution.rand)(x.parameters, x.gen);
}

///
auto randRange(T : SpecificDistribution!(ParametersType, Distribution, RGen),
                   ParametersType, Distribution, RGen)(T x)
    if (isSpecificDistribution!T)
{
    static import dstats.random;
    return dstats.random.randRange!(Distribution.rand)(x.parameters, x.gen);
}

///
@system unittest
{
    import std.random : Random, unpredictableSeed;
    import std.range : take;

    auto x = take(randRange!(Normal)(0, 2), 5);
    auto gen = Random(unpredictableSeed);
    auto y = take(randRange!(Normal)(0, 2, gen), 5);
}

///
@system unittest
{
    import std.range : take;

    auto stdNormal = specificDistribution!(Normal)(0, 1);

    auto x = take(randRange(stdNormal), 5);
}

@system unittest
{
    auto ber = randRange!Bernoulli(0.5);
    auto bet = randRange!Beta(1, 1);
    auto bin = randRange!Binomial(10, 0.5);
    auto cau = randRange!Cauchy(0, 2);
    auto chi = randRange!ChiSquare(2);
    auto exp = randRange!Exponential(1);
    auto fis = randRange!Fisher(1, 1);
    auto gam = randRange!Gamma(1, 1);
    auto geo = randRange!Geometric(0.5);
    auto hyp = randRange!Hypergeometric(2, 3, 4);
    auto lap = randRange!Laplace(0, 2);
    auto log = randRange!Logistic(2, 5);
    auto logn = randRange!LogNormal(0, 2);
    auto neg = randRange!NegBinom(10, 0.5);
    auto nor = randRange!Normal(0, 2);
    auto poi = randRange!Poisson(4);
    auto ray = randRange!Rayleigh(2);
    auto stu = randRange!StudentsT(10);
    auto wal = randRange!Wald(1, 1);
    auto wei = randRange!Weibull(1, 5);
}
