/*
TRAIT:
A traits class is class that is usually intended to be a meta-function associating 
types to other types or to constant values to provide a characterization of those types. 

POLICY:
A policy is normally meant to be a class that specifies what the behavior of another, 
generic class should be regarding certain operations that could be potentially 
realized in several different ways (and whose implementation is, therefore, left up 
to the policy class).


policies are passed as template parameters, while traits are global classes
*/

// todo: parallel threads, cuda
// grid trunctation 
// grid extension (for delta, gamma, theta)
// front fixing with BS control variate
// control variate in general

#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <numeric>
#include <array>
#include <chrono>
#include <type_traits>

/*
double normal_cdf(double x)
{
    double cdf;
    double poly;

    double xabs = std::abs(x);
    if (xabs > 37.0) 
        cdf = 0.0;
    else {  
        double exponential = exp( -xabs*xabs / 2.0);
        if (xabs < 7.07106781186547) { 
            poly = 3.52624965998911E-02 * xabs + 0.700383064443688;
            poly = poly * xabs + 6.37396220353165;
            poly = poly * xabs + 33.912866078383;
            poly = poly * xabs + 112.079291497871;
            poly = poly * xabs + 221.213596169931;
            poly = poly * xabs + 220.206867912376;
            cdf = exponential * poly;
            poly = 8.83883476483184E-02 * xabs + 1.75566716318264;
            poly = poly * xabs + 16.064177579207;
            poly = poly * xabs + 86.7807322029461;
            poly = poly * xabs + 296.564248779674;
            poly = poly * xabs + 637.333633378831;
            poly = poly * xabs + 793.826512519948;
            poly = poly * xabs + 440.413735824752;
            cdf = cdf / poly;
        } else {
            poly = xabs + 0.65;
            poly = xabs + 4 / poly;
            poly = xabs + 3 / poly;
            poly = xabs + 2 / poly;
            poly = xabs + 1 / poly;
            cdf = exponential / poly / 2.506628274631;
        }
    }
    if (x>0) cdf = 1.0 - cdf;
    return cdf;
}
*/

template<typename T>
void print(std::vector<T>& v) {
	for (auto i=v.begin(); i!=v.end(); ++i)
		std::cout << " " << *i;
	std::cout <<std::endl;	
}

/* ------------------------------------------------------------------------------------
 * traits and tags for derivatives
 * ------------------------------------------------------------------------------------ */ 
 // sub_type tags
 struct call_tag     {};
 struct put_tag      {};
 
 // exercise tags
 struct american_tag {};
 struct european_tag {};
 
template <class Derivative>
struct derivative_traits {
	typedef typename Derivative::sub_type sub_type;
	typedef typename Derivative::exercise_type exercise_type;
};

template <class Derivative>
double expiration_time(Derivative& d) {
	return d.T;
}

template <class Derivative>
double strike(Derivative& d) {
	return d.K;
}

// --------------------------------------------------------------------------------------
// vanilla_options
// --------------------------------------------------------------------------------------
template <typename SubType, typename ExerciseType>
struct vanilla_option {
	typedef SubType 		sub_type;
	typedef ExerciseType	exercise_type;
	double K;
	double T;
};


typedef vanilla_option<call_tag, european_tag> vanilla_european_call;
typedef vanilla_option<put_tag,  european_tag> vanilla_european_put;
typedef vanilla_option<call_tag, american_tag> vanilla_american_call;
typedef vanilla_option<put_tag,  american_tag> vanilla_american_put;

double payoff(vanilla_european_put& d, double S) {
	return std::max(d.K - S, 0.0);
}
double payoff(vanilla_american_put& d, double S) {
	return std::max(d.K - S, 0.0);
}
double payoff(vanilla_european_call& d, double S) {
	return std::max(S - d.K, 0.0);
}
double payoff(vanilla_american_call& d, double S) {
	return std::max(S - d.K, 0.0);
}

// ------------------------------------------------------------------------------------
// Stochastic processes
// ------------------------------------------------------------------------------------
struct Gbm {
	double sigma;
	double yield;
};

template <typename SDE>
double sigma(SDE& sde) { return sde.sigma; }

template <typename SDE>
double yield(SDE& sde) { return sde.yield; }

/* ------------------------------------------------------------------------------------
 * fill functions
 * ------------------------------------------------------------------------------------ */ 
template <typename Iterator>
void fill_exponential(Iterator b, Iterator e, double S, double lambda) {
	*b = S;
	std::generate(std::next(b), e, [&](){ S *= lambda; return S; } );
}

template <typename Iterator, typename Derivative>
void fill_exponential_payoff(Iterator b, Iterator e, double S, double lambda, Derivative& derivative) {
	*b = payoff(derivative,S);
	std::generate(std::next(b), e, [&](){ S *= lambda; return payoff(derivative,S); } );
}

/* ------------------------------------------------------------------------------------
 * generic binomial tree
 * ------------------------------------------------------------------------------------ */ 
template <typename Iterator, typename Derivative>
void binomial_tree_american_step(Iterator b, Iterator e, double wu, double wd, double S, double up_factor, Derivative& derivative) {
	for (Iterator i=b; i!=e; ++i) {
		*i = wd * (*i) + wu * (*(i+1));
		*i = std::max(*i, payoff(derivative,S) );
		S *= up_factor;
	}
}

template <typename Iterator, typename Derivative>
void binomial_tree_european_step(Iterator b, Iterator e, double wu, double wd, double /*S*/, double /*up_factor*/, Derivative& /*derivative*/) {
	for (Iterator i=b; i!=e; ++i) {
		*i = wd * (*i) + wu * (*(i+1));
		// No early exercise for European options
	}
}

template <typename Derivative>
double binomial_tree_price(double p, double u, double d, double discount_factor, double S0, Derivative& derivative, uint64_t n) {

	double wu = p * discount_factor;
	double wd = (1.0 - p) * discount_factor;
	double up_factor = u/d;	
	double S_bottom = S0 * std::pow(d,n);

	std::vector<double> V(n+1);
	fill_exponential_payoff(begin(V), end(V), S_bottom, up_factor, derivative);

	// Check if American or European exercise
	typedef typename Derivative::exercise_type exercise_type;
	constexpr bool is_american = std::is_same<exercise_type, american_tag>::value;

	for (uint64_t ti=n; ti>=1; --ti) {
		S_bottom /= d; // move left up
		if constexpr (is_american) {
			binomial_tree_american_step(&V[0], &V[ti], wu, wd, S_bottom, up_factor, derivative);
		} else {
			binomial_tree_european_step(&V[0], &V[ti], wu, wd, S_bottom, up_factor, derivative);
		}
	}
	return V[0];
}


// ------------------------------------------------------------------------------------
// Binomial tree types
// ------------------------------------------------------------------------------------

template <typename SDE, typename Derivative>
double CoxRossRubinsteinTreePrice(SDE& sde, double S0, Derivative& derivative, double r, uint32_t n) {
	double T = expiration_time(derivative);
	double dt = T / n;
	double u = exp( sigma(sde)*sqrt(dt));
	double d = exp(-sigma(sde)*sqrt(dt));
	double p = (exp(yield(sde) * dt) - d) / (u - d);
	double discount_factor = exp(-r*dt);
	return binomial_tree_price(p, u, d, discount_factor, S0, derivative, n);
}


template <typename SDE, typename Derivative>
double JarrowRuddTreePrice(SDE& sde, double S0, Derivative& derivative, double r, uint32_t n) {
	double T = expiration_time(derivative);
	double dt = T / n;
	double a = ( yield(sde) - 0.5*sigma(sde)*sigma(sde) ) * dt;
	double b = sigma(sde)*sqrt(dt);
	double u = exp(a + b);
	double d = exp(a - b);
	double p = 0.5;
	double discount_factor = exp(-r*dt);	
	return binomial_tree_price(p, u, d, discount_factor, S0, derivative, n);
}


template <typename SDE, typename Derivative>
double TianTreePrice(SDE& sde, double S0, Derivative& derivative, double r, uint32_t n) {
	double T = expiration_time(derivative);
	double dt = T / n;
	double vn = exp( sigma(sde) * sigma(sde) * dt);
	double u = 0.5 * exp(yield(sde) * dt) * vn * (vn + 1.0 + sqrt(vn*vn + 2.0*vn - 3.0));
	double d = 0.5 * exp(yield(sde) * dt) * vn * (vn + 1.0 - sqrt(vn*vn + 2.0*vn - 3.0));
	double p = (exp(yield(sde) * dt) - d) / (u - d);
	double discount_factor = exp(-r*dt);	
	return binomial_tree_price(p, u, d, discount_factor, S0, derivative, n);
}

/*
template <typename SDE>
struct TrigeorgisTree {
	binomial_params(double dt, SDE& sde, double& p, double& u, double* d) {
	}
};
*/
// http://www.wiwi.uni-bonn.de/sfb303/papers/1995/b/bonnsfb309.pdf
double LeisenReimer_f(double x, double n){
	return 0.5 + 0.5*(x>0 ? 1 : -1)*sqrt(1.0 - exp(-pow(x/(n+1.0/3.0+0.1/(n+1.0)),2)*(n+1.0/6.0)));
}

template <typename SDE, typename Derivative>
double LeisenReimerTree(SDE& sde, double S0, Derivative& derivative, double r, uint32_t n_) {
	uint32_t n = n_ + 1 - (n_%2); // force odd
	double T = expiration_time(derivative);
	double dt = T / n;
	double K = strike(derivative);
	double d1 = (log(S0/K) + (yield(sde) + 0.5*sigma(sde)*sigma(sde))*T) / (sigma(sde)*sqrt(T));
	double d2 = (log(S0/K) + (yield(sde) - 0.5*sigma(sde)*sigma(sde))*T) / (sigma(sde)*sqrt(T));
	double p = LeisenReimer_f(d2, n);
	double p_ = LeisenReimer_f(d1, n);
	double u = p_ / p * exp(yield(sde)*dt);
	//double d = (1.0 - p_) / (1.0 - p) * exp(yield(sde)*dt);
	double d = (exp(yield(sde)*dt) - p*u) / (1.0 - p);
	double discount_factor = exp(-r*dt);	
	return binomial_tree_price(p, u, d, discount_factor, S0, derivative, n);
}

// http://fbe.unimelb.edu.au/__data/assets/pdf_file/0010/806275/160.pdf
double ExtendedJoshi4_P4(double k,double dj)
{
	double alpha = dj/(sqrt(8.0));
	double alpha2 = alpha*alpha;
	double alpha3 = alpha*alpha2;
	double alpha5 = alpha3*alpha2;
	double alpha7 = alpha5*alpha2;
	double beta = -0.375*alpha - alpha3;
	double gamma = (5.0/6.0)*alpha5 + (13.0/12.0)*alpha3 +(25.0/128.0)*alpha;
	double delta = -0.1025*alpha - 0.9285*alpha3 - 1.43*alpha5 - 0.5*alpha7;
	double p =0.5;
	double rootk= sqrt(k);
	p+= alpha/rootk;
	p+= beta /(k*rootk);
	p+= gamma/(k*k*rootk);
	p+= delta/(k*k*k*rootk);
	return p;
}

template <typename SDE, typename Derivative>
double ExtendedJoshi4Tree(SDE& sde, double S0, Derivative& derivative, double r, uint32_t n_) {
	uint32_t n = n_ + 1 - (n_%2); // force odd
	double T = expiration_time(derivative);
	double dt = T / n;
	double K = strike(derivative);
	double d1 = (log(S0/K) + (yield(sde) + 0.5*sigma(sde)*sigma(sde))*T) / (sigma(sde)*sqrt(T));
	double d2 = (log(S0/K) + (yield(sde) - 0.5*sigma(sde)*sigma(sde))*T) / (sigma(sde)*sqrt(T));
	double p = ExtendedJoshi4_P4((n-1)/2, d2);
	double p_ = ExtendedJoshi4_P4((n-1)/2, d1);
	double u = p_ / p * exp(yield(sde)*dt);
	double d = (exp(yield(sde)*dt) - p*u) / (1.0 - p);
	double discount_factor = exp(-r*dt);	
	return binomial_tree_price(p, u, d, discount_factor, S0, derivative, n);
}



template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* ------------------------------------------------------------------------------------
 * main
 * ------------------------------------------------------------------------------------ */ 
int main() {

	{
		// Use same parameters as americanput.cpp for comparison
		double S0 = 100;
		double K = 100;
		double r = 0.05;
		double q = 0.05;
		double sigma = 0.25;
		double T = 1.0;
	std::cout.precision(16);	

		// Original test parameters (commented out for comparison with americanput.cpp)
		// S0 = 50;
		// K = 50;
		// r = 0.08;
		// q = 0.00;
		// sigma = 0.4;
		// T = 1.0;

		Gbm sde{sigma, r-q};
		vanilla_american_put american_derivative{K, T};
		vanilla_european_put european_derivative{K, T};


		uint64_t steps = 32;
		for (uint32_t i=0; i<10; ++i) {
			std::cout  << steps << " ";
			
			// European put price
			{
				auto t1 = std::chrono::high_resolution_clock::now();
				double  price = CoxRossRubinsteinTreePrice(sde, S0, european_derivative, r, steps);
				auto t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
				std::cout << "EU:" << price << "(" << fp_ms.count() << "ms) ";
			}
			
			// American put prices
			{
				auto t1 = std::chrono::high_resolution_clock::now();
				double  price = CoxRossRubinsteinTreePrice(sde, S0, american_derivative, r, steps);
				auto t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
				std::cout << "AM_CRR:" << price << "(" << fp_ms.count() << "ms) ";
			}

			{
				auto t1 = std::chrono::high_resolution_clock::now();
				double  price = TianTreePrice(sde, S0, american_derivative, r, steps);
				auto t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
				std::cout << "AM_Tian:" << price << "(" << fp_ms.count() << "ms)";
			}
 
 			std::cout  << std::endl;
 			
			steps *= 2;
		}	
		return 0;
	}
	
	std::cout.precision(std::numeric_limits<long double>::max_digits10);	
}