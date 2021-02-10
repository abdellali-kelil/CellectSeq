//---------------------------------------------------------------------------
//
// boost library
#include <boost/math/distributions/students_t.hpp>

// ---------------------------------------------------------------------------
// Effect size based on difference between means (Cohen's index)
double Effect(double avr1, double avr2, double var1, double var2, double nbr1, double nbr2)
{
	/*
		Effect size   |   d    |  Reference
		==============|========|===================
		Very small    |  0.01  |  Sawilowsky, 2009
		Small         |  0.20  |  Cohen, 1988
		Medium        |  0.50  |  Cohen, 1988
		Large         |  0.80  |  Cohen, 1988
		Very large    |  1.20  |  Sawilowsky, 2009
		Huge          |  2.00  |  Sawilowsky, 2009
	*/

	// Calculate Cohen's d
	return (avr1 - avr2) / sqrt( ((nbr1 - 1) * var1 + (nbr2 - 1) * var2) / (nbr1 + nbr2 - 2));
}

// ---------------------------------------------------------------------------
// Calculate the p-value of the t-test for comparing the means of two samples
// See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
// ---------------------------------------------------------------------------
double pvalue(double avr1, double avr2, double sdv1, double sdv2, double nbr1, double nbr2)
{
	/*
		avr1 = Sample 1 Mean
		sdv1 = Sample 1 SD
		nbr1 = Sample 1 Size
		avr2 = Sample 2 Mean
		sdv2 = Sample 2 SD
		nbr2 = Sample 2 Size
	*/

	// check nbr1 and nbr2
	if (nbr1 <= 1 || nbr2 <= 1) return 1.0f;

	// check sdv1 and sdv2
	if (sdv1 == 0 && sdv2 == 0) return 1.0f;

	// check avr1 and avr2
	if (avr1 == avr2) return 1.0f;

	// degrees of freedom
	double df = sdv1 * sdv1 / nbr1 + sdv2 * sdv2 / nbr2; df *= df;
	double t1 = sdv1 * sdv1 / nbr1; t1 *= t1; t1 /= (nbr1 - 1);
	double t2 = sdv2 * sdv2 / nbr2; t2 *= t2; t2 /= (nbr2 - 1);

	// normalize
	df /= (t1 + t2);

	// t-score
	double ts = (avr1 - avr2) / sqrt(sdv1 * sdv1 / nbr1 + sdv2 * sdv2 / nbr2);

	// p-value
	try
	{
		// Prepare distribution
		boost::math::students_t dis(df);

		// Calculate pvalue for specificity
		return boost::math::cdf(complement(dis, ts));

		// Check error
		printf_s("FATAL EEROR !!!\n"), exit(EXIT_FAILURE);
	}

	catch(std::exception& e)
	{
		// p-value not calculated
		printf_s("\nError p-value %s\n", e.what()); exit(EXIT_FAILURE);
	}
}
