#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include "Poisson.h"
using std::exp;
using std::floor;
using std::max;
using std::ceil;
using std::sqrt;
int FoxFindRTP(double lambda, double epsilon, int m, double &k_r)
{
	double a_lambda;
	double k;
	double d_k_lambda;
	int R;

	/* Fox-Glynn Page 443, Column 1 Rule 2 */
	if (lambda < 400.0)
		lambda = 400.0;

	a_lambda = (1.0 + 1.0 / lambda) * exp(1.0 / 16.0) * sqrt(2.0);

	for (k = 3.0;; k += 1.0)
	{
		d_k_lambda = 1.0 / (1.0 - exp((-2.0 / 9.0) * (k * sqrt(2.0 * lambda) + 1.5)));
		if (a_lambda * d_k_lambda * exp(-k * k / 2.0) / (k * sqrt(2.0 * M_PI)) < epsilon / 2.0)
		{
			R = (int)ceil(m + k * sqrt(2.0 * lambda) + 1.5);
			k_r = k;
			break;
		}
	}
	return R;

}

static double TAU = 1e-60;
static double OMEGA = 1e60;
void fox_find_trunc_point(double lambda, int &Ltp, int &Rtp, double epsilon, bool &overflowed)
{
	/************************************************************************
	finds the left truncation point and calls the routine foc_find_r to
	find the right truncation point

	lambda =  qt, basically the Poisson variable for Poisson probability
	computation. It is an input parameter.

	epsilon = error tolerance

	Ltp, Rtp = output parameters, carry the value of LTP and RTP

	*/

	double k;
	double k_r;
	int L, R;
	double b_lambda;
	double c_m;
	double k_hat;
	double k_tilde;
	double bound = 0.0;
	int m = floor(lambda);
	
	overflowed = false;
	Rtp = R = FoxFindRTP(lambda, epsilon, m, k_r);

	if (lambda < 25.0)
	{
		Ltp = 0;
		return;
	}
	b_lambda = (1.0 + 1.0 / lambda) * exp(1.0 / (8.0 * lambda));

	for (k = 1.0;; k += 1.0)
	{
		if (b_lambda * exp(-k * k / 2.0) / (k * sqrt(2.0 * M_PI)) < epsilon / 2.0)
		{
			break;
		}
	}
	L = (int)floor(m - k * sqrt(lambda) - 1.5);
	if (L < 0)
		L = 0;
	Ltp = L;
	/*
	*   The following code checks the lower bounds on the
	*   poisson probabilities to see if underflow can occur.
	*/
	c_m = (1.0 / sqrt(2.0 * M_PI * m)) * exp((m - lambda - 1) / (12.0 * m));

	if (lambda >= 400.0)
	{
		k_hat = k_r * sqrt(2.0) + 3.0 / (2.0 * sqrt(lambda));
		bound = c_m * exp(-(k_hat + 1) * (k_hat + 1) / 2.0);
		if (OMEGA * bound / (1.0e10 * (R - L)) < TAU)
		{
			overflowed = true;
		}
	} /* if (lambda >= 400.0) */

	k_tilde = k + 3.0 / (2.0 * sqrt(lambda));

	if (k_tilde > 0 && k_tilde <= sqrt(lambda) / 2.0)
	{
		bound = c_m * exp(-k_tilde * k_tilde / 2.0 -
			-k_tilde * k_tilde * k_tilde / (3.0 * sqrt(lambda)));
	}
	else if (k_tilde <= sqrt((double)(m + 1)))
	{
		bound = max(c_m *pow(1.0 - k_tilde / sqrt((double)(m + 1)),
			k_tilde * sqrt((double)(m + 1))), exp(-lambda));
	}
	if (bound != 0.0 && OMEGA * bound / (1.0e10 * (R - L)) < TAU)
	{
		overflowed = true;
	}
}

double fox_weight_sum(const std::vector<double>& weight)
{
	/* Compute W, to attenuate roundoff, small terms are added first.
	W = w[L] + .... + w[R]  */
	double W = 0;
	int s = 0;
	int t = weight.size() - 1;

	while (s < t)
	{
		if (weight[s] <= weight[t])
		{
			W += weight[s];
			s++;
		}
		else
		{
			W += weight[t];
			t--;
		}
	}
	W += weight[s];
	return W;
}

std::vector<double> fox_term_weight(double lambda, int L, int R, double epsilon, bool &overflowed)
{

	/************************************************************************
	lambda = qt , input parameter
	L, R = LTP and RTP, input parameters

	W = a large weight assigned to avoid underflow, output parameter.

	w = array of weights. w[i] = weight for ith term where L <= i <= R

	General Comments

	epsilon = error tolerance, greater than 10^(-10).
	TAU = Underflow threshold
	OMEGA = overflow threshold
	w(i) = p(i)/W  where p(i) is the Poisson probability.
	W = Total wt.
	for lambda >= 400, mass left of L <= epsilon/2
	mass right of L <= epsilon/2
	for 0 < lambda < 400, mass left of L <= epsilon/2
	mass right of L <= epsilon/2 + 6*10^12 *TAU/OMEGA
	oflow_flag = 1 if no underflow can occur while computing the weights
	oflow_flag = 0 indicates that wts are not computed due to potential underflow.
	w(i)/W may underflow even if F = 1
	*/

	int j, s, t;
	double qty;
	std::vector<double> weight(R - L + 1);
	int m = floor(lambda);

	weight[m - L] = OMEGA / (1.0e10 * (R - L));

	/* down part from Fox-Glynn Page 442 */
	for (j = m; j > L; j--)
		weight[j - 1 - L] = (j / lambda) * weight[j - L];

	if (lambda < 400.0)
	{
		if (R > 600)
		{
			overflowed = true;
		}
		for (j = m; j < R; j++)
		{
			qty = lambda / (j + 1);
			if (weight[j - L] > TAU / qty)
				weight[j + 1 - L] = qty * weight[j - L];
			else
				R = j;
		}
	}
	else
	{
		for (j = m; j < R; j++)
			weight[j + 1 - L] = weight[j - L] * (lambda / (j + 1));
	}
	
	return weight;
}

int find_cum_rtp(double lambda, double precision, double time)
{
	double log_qt = log(lambda);
	double tmpti = -lambda;
	double p_term = exp(tmpti);
	double sumqt = p_term;
	double right = 1;

	while (p_term < TAU)
	{
		tmpti += (log_qt - log(right));
		p_term = exp(tmpti);
		sumqt += p_term;
		right += 1;
	} /* end while */

	while ((precision < time * (1.0 - sumqt)) && (p_term > TAU))
	{
		tmpti += (log_qt - log(right));
		p_term = exp(tmpti);
		sumqt += p_term;
		right += 1;
	} /* end while */

	return right;
}
