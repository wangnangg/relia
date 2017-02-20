#pragma once
#include <vector>


double fox_weight_sum(const std::vector<double>& weight);
void fox_find_trunc_point(double lambda, int &Ltp, int &Rtp, double epsilon, bool &overflowed);
std::vector<double> fox_term_weight(double lambda, int L, int R, double epsilon, bool &overflowed);

int find_cum_rtp(double lambda, double precision, double time);
