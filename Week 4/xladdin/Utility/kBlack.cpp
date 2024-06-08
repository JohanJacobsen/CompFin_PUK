#include "kBlack.h"
#include "kSolver.h"
#include "kFd1d.h"

class kBlackObj : public kSolverObjective
{
public:

	//	constructor
	kBlackObj(
		double	expiry,
		double	strike,
		double	price,
		double	forward)
		: kSolverObjective(),
		myExpiry(expiry),
		myStrike(strike),
		myPrice(price),
		myForward(forward)
	{}

	//	value
	virtual double	value(double x)
	{
		double res = kBlack::call(myExpiry, myStrike, myForward, x) - myPrice;

		//	done
		return res;
	}

	//	deriv
	virtual double	deriv(double x)
	{
		double res = kBlack::vega(myExpiry, myStrike, myForward, x);

		//	done
		return res;
	}

	//	private parts
private:

	//	expiry
	double	myExpiry;
	double	myStrike;
	double	myPrice;
	double	myForward;

};

// implied vol
double
kBlack::implied(
	double	expiry,
	double	strike,
	double	price,
	double	forward)
{
	//	calc intrinsic
	double intrinc = max(forward - strike, 0.0);
	if (price <= intrinc) return 0.0;

	//	objective
	kBlackObj obj(expiry, strike, price, forward);

	//	start guess
	double volatility = 0.1;
	int    numIter = 10;
	double epsilon = (price - intrinc) * kConstants::epsilon();

	//	solve
	kSolver::newtonRapson(obj, volatility, numIter, epsilon, nullptr);

	//	bound
	volatility = max(0.0, volatility);

	//	done
	return volatility;
}

bool
kBlack::fdRunner(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	const int			boundary, //	don't calc boundary (0), calc boundary (1)
	const double		barrier,    //  down-and-out barrier
	const double		bt, // down-and-out type simple (0), simple + barrier point (1), or advanced (2)
	double& res0,
	kVector<double>& s,
	kVector<double>& res,
	kVector<double>& eeb, // output vector for Early exercise boundary
	string& error)
{
	//	helps
	int h, i, p;

	//	construct s axis
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double xl = -numStd * std;
	double xu = numStd * std;
	int    nums = 2 * (numS / 2);
	double dx = (xu - xl) / max(1, nums);
	if (nums <= 0 || xl == xu)
	{
		nums = 1;
	}
	else
	{
		++nums;
	}

	s.resize(nums);
	double ds = exp(dx);
	s(0) = s0 * exp(xl);
	for (i = 1; i < nums; ++i)
	{
		s(i) = s(i - 1) * ds;
	}
	if (bt > 0)
	{
		if (barrier > 1) {
			// Find the index where the barrier should be inserted
			int barrierIndex = -1;
			for (i = 0; i < nums; ++i)
			{
				if (s(i) > barrier) {
					barrierIndex = i;
					break;
				}
			}
			if (barrierIndex != -1) {
				// Shift the elements to the right to make space for the barrier
				s.resize(nums + 1);
				for (i = nums; i > barrierIndex; --i) {
					s(i) = s(i - 1);
				}
				// Insert the barrier at the appropriate index
				s(barrierIndex) = barrier;
				nums++; // Increment the number of points
			}
		}
	}



	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, s, false);

	//	set terminal result
	double sl, su;
	res.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		if (s(i) <= barrier) {
			res(i) = 0.0;
		}
		else {
			if (smooth == 0 || i == 0 || i == nums - 1)
			{
				if (dig) res(i) = 0.5 * (kInlines::sign(s(i) - strike) + 1.0);
				else    res(i) = max(0.0, s(i) - strike);
			}
			else
			{
				sl = 0.5 * (s(i - 1) + s(i));
				su = 0.5 * (s(i) + s(i + 1));
				if (dig) res(i) = kFiniteDifference::smoothDigital(sl, su, strike);
				else	res(i) = kFiniteDifference::smoothCall(sl, su, strike);
			}

			if (pc < 0)
			{
				if (dig) res(i) = 1.0 - res(i);
				else    res(i) -= (s(i) - strike);
			}
		}
	}
	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	//  initialize eeb vector
	eeb.resize(numT + 1);
	eeb(numt) = strike;

	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		if (bt > 0)
		{
			for (i = 0; i < nums; ++i)
			{
				if (s(i) <= barrier) {
					fd.r()(i) = r;
					fd.mu()(i) = 0.0;
					fd.var()(i) = 0.0;
				}
				else {
					fd.r()(i) = r;
					fd.mu()(i) = mu * s(i);
					fd.var()(i) = kInlines::sqr(sigma * s(i));
				}
			}
		}
		else
		{
			for (i = 0; i < nums; ++i)
			{
				fd.r()(i) = r;
				fd.mu()(i) = mu * s(i);
				fd.var()(i) = kInlines::sqr(sigma * s(i));
			}
		}


		//	roll
		fd.res()(0) = res;
		for (h = numt - 1; h >= 0; --h)
		{
			// Apply the barrier condition after each backward step
			fd.rollBwd(dt, update || h == (numt - 1), theta, wind, fd.res());
			for (i = 0; i < nums; ++i)
			{
				if (s(i) <= barrier) {
					fd.res()(0)(i) = 0.0; // Set to zero if barrier is crossed
				}
			}
			//Early exercise boundary
			double exercise = s(0);
			if (ea > 0)
			{
				for (i = 0; i < nums; ++i)
				{
					fd.res()(0)(i) = max(res(i), fd.res()(0)(i));
					if (boundary > 0)
					{
						if (pc > 0) {
							if (res(i) != fd.res()(0)(i))
							{
								exercise = s(i + 1);
							}
						}
						else {
							if (res(i) == fd.res()(0)(i)) {
								exercise = s(i);
							}
						}
					}
				}
				eeb(h) = exercise;
			}
		}
	}

	//	set result
	res = fd.res()(0);
	res0 = fd.res()(0)(nums / 2);

	//	done
	return true;
}

bool
kBlack::fdfwdRunner(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double        max_expiry,
	kVector<double>& strikes,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	//double&			res0,
	kVector<double>& s,
	kMatrix<double>& resM,
	string& error)
{
	int numK = strikes.size();
	resM.resize(numT + 1, numK);

	double t = max(0.0, max_expiry);
	double std = sigma * sqrt(t);
	double sl = log(s0) - numStd * std;
	double su = log(s0) + numStd * std;

	int    nums = 2 * (numS / 2);
	double ds = (su - sl) / max(1, nums);
	if (nums <= 0 || sl == su)
	{
		nums = 1;
	}
	else
	{
		++nums;
	}
	s.resize(nums);

	s(0) = sl;
	for (int i = 1; i < nums; ++i)
	{
		s(i) = s(i - 1) + ds;
	}

	for (int i = 0; i < nums; ++i) {
		s(i) = exp(s(i));
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, s, false);

	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	for (int k = 0; k < numK; k++) {

		double strike = strikes[k];
		kVector<double> payoff(nums);
		kVector<double> res(nums);

		for (int i = 0; i < nums; ++i) {
			res(i) = 0;
			payoff(i) = max(0.0, s(i) - strike);
		}
		res(nums / 2) = 1;
		resM(0, k) = payoff(nums / 2) * res(nums / 2);

		// Forward rolling 
		for (int h = 0; h < numt; ++h) {
			// Set parameters 
			for (int i = 0; i < nums; ++i) {
				fd.r()(i) = r;
				fd.mu()(i) = mu * s(i);
				fd.var()(i) = sigma * sigma * s(i) * s(i);
			}
			// Perform the forward roll 
			fd.res()(0) = res;
			fd.rollFwd(dt, update || h == 0, theta, wind, fd.res());

			res = fd.res()(0);

			for (int i = 0; i < nums; ++i) {
				resM(h + 1, k) += payoff(i) * res(i);
			}
		}
	}
	return true;
}

bool
kBlack::fdfwdRunner_call(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double        max_expiry,
	//kVector<double>&	strikes,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	//double&			res0,
	kVector<double>& s,
	kMatrix<double>& resM,
	string& error)
{
	int numK = 2 * (numS / 2);
	resM.resize(numT + 1, numK + 1);

	double t = max(0.0, max_expiry);
	double std = sigma * sqrt(t);

	double sl = log(s0) - numStd * std;
	double su = log(s0) + numStd * std;

	double ds = (su - sl) / max(1, numK);
	if (numK <= 0 || sl == su)
	{
		numK = 1;
	}
	else
	{
		++numK;
	}
	s.resize(numK);

	s(0) = sl;
	for (int i = 1; i < numK; ++i)
	{
		s(i) = s(i - 1) + ds;
	}

	for (int i = 0; i < numK; ++i) {
		s(i) = exp(s(i));
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, s, false);

	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	kVector<double> res(numK);

	for (int k = 0; k < numK; k++) {
		res(k) = max(0.0, s0 - s(k));

		resM(0, k) = res(k);
	}
	// Forward rolling 
	for (int h = 0; h < numt; ++h) {
		// Set parameters 
		for (int i = 0; i < numK; ++i) {
			fd.r()(i) = r;
			fd.mu()(i) = mu * s(i);
			fd.var()(i) = sigma * sigma * s(i) * s(i);
		}
		// Perform the call roll 
		fd.res()(0) = res;
		fd.rollCall(dt, update || h == 0, theta, wind, fd.res());

		res = fd.res()(0);

		for (int i = 0; i < numK; ++i) {
			resM(h + 1, i) = res(i);
		}
	}
	return true;
}

//	fd runner
bool
kBlack::fdRunnerTry(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	double& res0,
	kVector<double>& s,
	kVector<double>& res,
	string& error)
{
	//	helps
	int h, i, p;

	//	construct s axis
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double sl = log(s0) - numStd * std;
	double su = log(s0) + numStd * std;
	int    nums = 2 * (numS / 2);
	double ds = (su - sl) / max(1, nums);
	if (nums <= 0 || sl == su)
	{
		nums = 1;
	}
	else
	{
		++nums;
	}
	s.resize(nums);

	s(0) = sl;
	for (i = 1; i < nums; ++i)
	{
		s(i) = s(i - 1) + ds;
	}

	kVector<double> payoff(nums);

	for (i = 0; i < nums; ++i) {
		s(i) = exp(s(i));
		payoff(i) = max(0.0, s(i) - strike);
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, s, false);

	//	set terminal result
	res.resize(nums);

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);


	//	set parameters
	for (i = 0; i < nums; ++i)
	{
		fd.r()(i) = r;
		fd.mu()(i) = mu * s(i);
		fd.var()(i) = sigma * sigma * s(i) * s(i);
	}

	//initial condition 
	for (int i = 0; i < nums; ++i) {
		res(i) = 0;
	}

	res(nums / 2) = 1;
	//	roll
	fd.res()(0) = res;

	for (int h = 0; h < numt; ++h) {
		fd.rollFwd(dt, update || h == 0, theta, wind, fd.res());
	}

	res = fd.res()(0);
	//res0 = fd.res()(0)(nums / 2);

	res0 = 0.0;
	for (int i = 0; i < nums; ++i) {
		res0 += payoff(i) * res(i);
	}

	return true;
}

bool
kBlack::fdfwdRunner_new(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	kVector<double>& expiries,
	kVector<double>& strikes,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	kVector<double>& s,
	kMatrix<double>& resM,
	string& error)
{

	int numK = strikes.size();
	int numE = expiries.size();

	resM.resize(numE + 1, numK); //Muligvis uden +1

	for (int e = 0; e < numE; ++e) {

		double t = max(0.0, expiries[e]);
		double std = sigma * sqrt(t);

		double sl = log(s0) - numStd * std;
		double su = log(s0) + numStd * std;

		int    nums = 2 * (numS / 2);
		double ds = (su - sl) / max(1, nums);
		if (nums <= 0 || sl == su)
		{
			nums = 1;
		}
		else
		{
			++nums;
		}
		s.resize(nums);

		s(0) = sl;
		for (int i = 1; i < nums; ++i)
		{
			s(i) = s(i - 1) + ds;
		}

		for (int i = 0; i < nums; ++i) {
			s(i) = exp(s(i));
		}

		//	construct fd grid
		kFd1d<double> fd;
		fd.init(1, s, false);

		int    numt = max(0, numT);
		double dt = t / max(1, numt);

		for (int k = 0; k < numK; k++) {

			double strike = strikes[k];

			kVector<double> payoff(nums);
			kVector<double> res(nums);

			for (int i = 0; i < nums; ++i) {
				res(i) = 0;
				payoff(i) = max(0.0, s(i) - strike);
			}

			res(nums / 2) = 1;
			resM(0, k) = payoff(nums / 2) * res(nums / 2);

			// Forward rolling 
			for (int h = 0; h < numt; ++h) {
				// Set parameters 
				for (int i = 0; i < nums; ++i) {
					fd.r()(i) = r;
					fd.mu()(i) = mu * s(i);
					fd.var()(i) = sigma * sigma * s(i) * s(i);
				}

				// Perform the forward roll 
				fd.res()(0) = res;
				fd.rollFwd(dt, update || h == 0, theta, wind, fd.res());

				res = fd.res()(0);
			}

			for (int i = 0; i < nums; ++i) {
				resM(e + 1, k) += payoff(i) * res(i);
			}
		}
	}
	return true;
}

