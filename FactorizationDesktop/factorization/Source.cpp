#pragma once
#include "Header.h"
#include "..\FactorizationDesktop.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "AtkinSieve.cpp"
#include <map>
#include <set>
#include <mpirxx.h>
#include <cstdint>
#include <mpfr.h>
#include <numeric>
#pragma comment (lib, "mpir.lib")
#pragma comment (lib, "mpirxx.lib")
using namespace std;
namespace fact {
	int mil_rab(mpz_class n) // 1-prostoe
	{

		mpz_class k, s = 0, buf = 0;;
		mpz_class t, tt, a, x, i;
		k = 10; //к-во итераций

		if (n <= 3)
		{
			return 1;
		}

		t = n - 1;

		//while (t % 2 == 0) {  //представление n-1 = (2^s)t, где t нечетно
		while (mpz_even_p(t.get_mpz_t())) {
			t = t / 2;
			s++;
		}

		for (mpz_class i = 0; i < k; i++)
		{
			a = rand();
			a = a % (n - 2) + 2; //случайное число из [2,n-1]

			//x = powmod(a, t, n);
			mpz_powm(x.get_mpz_t(), a.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
			if ((x == 1) || (x == n - 1)) continue;
			for (mpz_class j = 1; j < s; j++)
			{
				//x = pow(x, 2) % n;
				//x = powmod(x, 2, n);
				mpz_powm(x.get_mpz_t(), x.get_mpz_t(), mpz_class(2).get_mpz_t(), n.get_mpz_t());
				if (x == 1)
				{
					return 0;
				}
				if (x == n - 1) break;
			}
			if (x != n - 1)
			{
				return 0;
			}
		}
		return 1;
	}

	int Miller(mpz_class n)
	{
		mpz_class f = sqrt(sqrt(sqrt(n))) * 2;
		mpz_class a = 2;
		mpz_class buf;
		bool fl = false;
		mpz_class v, k = 0;
		for (a; a <= f; a = a + 1)
		{
			if (n%a == 0) return 0;
			mpz_class temp;
			mpz_powm(temp.get_mpz_t(), a.get_mpz_t(), mpz_class(n - 1).get_mpz_t(), n.get_mpz_t());
			if (temp != 1) return 0;

			buf = 2;
			k = 0;
			while (buf < n)
			{
				if ((n - 1) % buf == 0) v = k;
				k++;
				buf = buf * 2;
			}
			for (k = 1; k <= v; k++)
			{
				//buf = powmod(a, (n - 1) / pow((int)2, k), n);
				mpz_class power;
				mpz_pow_ui(power.get_mpz_t(), mpz_class(2).get_mpz_t(), k.get_ui());
				power = (n - 1) / power;
				mpz_powm(buf.get_mpz_t(), a.get_mpz_t(), power.get_mpz_t(), n.get_mpz_t());
				if (buf == 0) buf = n - 1;
				else buf = buf - 1;

				//buf = Nod(buf, n);
				mpz_gcd(buf.get_mpz_t(), buf.get_mpz_t(), n.get_mpz_t());
				if ((buf > 1) && (buf < n))
				{
					fl = true;
					break;
				}
			}
			if (fl == true) return 0;
		}
		return 1;
	}
	/*
	mpz_class LegendreSymbol(mpz_class a, mpz_class p) {
		if (a >= p)
			a %= p;
		mpz_class res;// = powmod(a, (p - 1) / 2, p);
		mpz_powm(res.get_mpz_t(), a.get_mpz_t(), mpz_class((p-1)/2).get_mpz_t(), p.get_mpz_t());
		return res > 1 ? -1 : res;
	}
	*/
	mpz_class SqrtMod(mpz_class a, mpz_class p) {
		if (a >= p)
			a %= p;
		mpz_class n = 1;
		//while (LegendreSymbol(n, p) >= 0) {
		while (mpz_legendre(n.get_mpz_t(), p.get_mpz_t()) >= 0) {
			++n;
		}
		mpz_class alpha = 0;
		mpz_class s = p - 1;
		//while (s % 2 == 0) {
		while (mpz_even_p(s.get_mpz_t())) {
			++alpha;
			s /= 2;
		}
		mpz_class b;// = powmod(n, s, p);
		mpz_powm(b.get_mpz_t(), n.get_mpz_t(), s.get_mpz_t(), p.get_mpz_t());
		mpz_class r;// = powmod(a, (s + 1) / 2, p);
		mpz_powm(r.get_mpz_t(), a.get_mpz_t(), mpz_class((s + 1) / 2).get_mpz_t(), p.get_mpz_t());
		mpz_class rCalc;// = powmod(a, s, p);
		mpz_powm(rCalc.get_mpz_t(), a.get_mpz_t(), s.get_mpz_t(), p.get_mpz_t());
		mpz_class j;
		mpz_class power;
		mpz_pow_ui(power.get_mpz_t(), mpz_class(2).get_mpz_t(), mpz_class(alpha - 2).get_ui());
		mpz_class check;// = powmod(rCalc, pow(2, alpha - 2), p);
		mpz_powm(check.get_mpz_t(), rCalc.get_mpz_t(), power.get_mpz_t(), p.get_mpz_t());
		if (check > 1)
			check -= p;
		if (check == 1) {
			j = 0;
		}
		else
			if (check == -1) {
				j = 1;
			}
			else
				cout << "check error" << endl;
		for (mpz_class i = 1; i < alpha - 1; ++i) {
			//check = powmod(powmod(powmod(b, j, p), 2, p) * rCalc, pow(2, alpha - i - 2), p);
			mpz_class base, exp;
			mpz_pow_ui(exp.get_mpz_t(), mpz_class(2).get_mpz_t(), mpz_class(alpha - i - 2).get_ui());
			mpz_powm(base.get_mpz_t(), b.get_mpz_t(), j.get_mpz_t(), p.get_mpz_t());
			mpz_powm(base.get_mpz_t(), base.get_mpz_t(), mpz_class(2).get_mpz_t(), p.get_mpz_t());
			base = base * rCalc;
			mpz_powm(check.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
			if (check > 1)
				check -= p;
			if (check == 1) {
				//j = 0;
			}
			else
				if (check == -1) {
					mpz_class addition;
					mpz_pow_ui(addition.get_mpz_t(), mpz_class(2).get_mpz_t(), i.get_ui());
					//j += pow(2, i);
					j += addition;
				}
				else
					std::cout << "check error" << endl;
		}
		mpz_class result;
		mpz_powm(result.get_mpz_t(), b.get_mpz_t(), j.get_mpz_t(), p.get_mpz_t());
		//return (powmod(b, j, p) * r) % p;
		return (result * r) % p;
	}
	/*
	mpz_class gcd(mpz_class a, mpz_class b) {
		while (b != 0) {
			a %= b;
			swap(a, b);
		}
		return a;
	}

	mpz_class gcdex(mpz_class a, mpz_class b, mpz_class & x, mpz_class & y) {//–асширенный алгоритм евклида emaxx
		if (a == 0) {
			x = 0; y = 1;
			return b;
		}
		mpz_class x1, y1;
		mpz_class d = gcdex(b%a, a, x1, y1);
		x = y1 - (b / a) * x1;
		y = x1;
		return d;
	}

	mpz_class ReverseMod(mpz_class a, mpz_class m) {//обратное по модулю число
		mpz_class x = 0, y;
		mpz_class g = gcdex(a, m, x, y);
		if (g != 1)
			std::cout << "ReverseMod no solution\n";
		else {
			x = (x % m + m) % m;
		}
		return x;
	}
	*/
	mpz_class sqrtModPowerLiftingPlusOne(mpz_class root, mpz_class quadResidue, mpz_class mod, uint64_t modBase) {
		mpz_class newmod = mod * modBase;
		//mpz_class yk = ReverseMod(2 * root, mod)% newmod;	
		mpz_class yk;
		mpz_class r = 2 * root;
		if (mpz_invert(yk.get_mpz_t(), r.get_mpz_t(), mod.get_mpz_t()) == 0)
			std::cout << "ReverseMod no solution\n";
		yk %= newmod;
		quadResidue %= newmod;
		root = root - (((root*root) % newmod - quadResidue)*yk) % newmod;

		if (root < 0) {
			root = -root;
		}
		root %= newmod;
		//root = min(root, -(root - newmod));
		return root;
	}

	bool chekNuberInRange(mpz_class lowerBound, mpz_class upperBound, mpz_class t, mpz_class p, mpz_class A) {// check if there is range lower <=T<= upper that T = t1 or t2 ( mod p)
		if (A >= p) {
			return true;
		}
		mpz_class l = lowerBound % p;
		mpz_class u = upperBound % p;
		if (l >= u) {
			if (t >= u || t <= l)
				return true;
		}
		else {
			if (t >= l && t <= u)
				return true;
		}
		return false;
	}

	uint64_t Step4(mpz_class n, uint64_t p, mpz_class A, mpz_class& res1, mpz_class& res2) {
		mpz_class lowerBound = sqrt(n) + 1;
		mpz_class upperBound = sqrt(n) + A;

		mpz_class residue = n;
		uint64_t modBase = p;
		mpz_class currentMod = p;

		mpz_class root = SqrtMod(residue, currentMod).get_ui();
		mpz_class root2 = -(root - modBase);
		res1 = root;
		res2 = root2;
		uint64_t b = 1;
		while (!chekNuberInRange(lowerBound, upperBound, root, currentMod, A)
			&& !chekNuberInRange(lowerBound, upperBound, root2, currentMod, A)
			&& (res1 <= upperBound || res2 <= upperBound)) {
			root = sqrtModPowerLiftingPlusOne(mpz_class(root), residue, currentMod, modBase).get_ui();
			currentMod *= modBase;
			root2 = -(root - currentMod.get_ui());
			++b;
			res1 = root;
			res2 = root2;
		}
		if (res1 > upperBound && res2 > upperBound)
			return 0;
		do {
			res1 = root;
			res2 = root2;
			root = sqrtModPowerLiftingPlusOne(mpz_class(root), residue, currentMod, modBase).get_ui();
			currentMod *= modBase;
			root2 = -(root - currentMod.get_ui());
			++b;
		} while (chekNuberInRange(lowerBound, upperBound, root, currentMod, A)
			|| chekNuberInRange(lowerBound, upperBound, root2, currentMod, A));
		return b - 1;
	}

	mpz_class liftRootEvenModulo(mpz_class root, mpz_class a, mpz_class p) {
		//mpz_class i = ((root * root - a) / p) % 2;
		//return root + i * (p / 2);
		mpz_class buf((root * root - a) / p);
		if (mpz_odd_p(buf.get_mpz_t()))
			return root + p / 2;
		else
			return root;
	}

	class GaussianEliminator {
		vector<mpz_class> listingT;
		vector<uint64_t> primes;
		uint64_t lastBasisPos = 0;
		uint64_t primeBaseSize;
		vector<uint64_t> fillOrder;
		set<uint64_t> remainingPositions;
		vector<vector<uint64_t>> operationsListVector;
		vector<vector<uint64_t>> exponentBaseMatrix;
	public:
		GaussianEliminator(vector<uint64_t>& factorBase) {
			primeBaseSize = factorBase.size();
			listingT = vector<mpz_class>(primeBaseSize);
			primes = factorBase;
			fillOrder = vector<uint64_t>(primeBaseSize);
			exponentBaseMatrix = vector<vector<uint64_t>>(primeBaseSize);
			remainingPositions = set<uint64_t>();
			operationsListVector = vector<vector<uint64_t>>(primeBaseSize, vector<uint64_t>());
			for (uint64_t i = 0; i < primeBaseSize; ++i)
				remainingPositions.insert(i);
		}

		mpz_class addVector(mpz_class& n, mpz_class& T, vector<uint64_t>& exponentVector) {

			//eliminate vector
			vector<uint64_t> eliminatedVector(primeBaseSize);
			for (uint64_t i = 0; i < primeBaseSize; ++i)
				eliminatedVector[i] = exponentVector[i] % 2;

			for (uint64_t i = 0; i < lastBasisPos; ++i) {
				uint64_t pos = fillOrder[i];
				if (eliminatedVector[pos] == 0)
					continue;
				vector<uint64_t>& operations = operationsListVector[pos];
				uint64_t operationsCount = operations.size();
				for (uint64_t j = 0; j < operationsCount; ++j) {
					eliminatedVector[operations[j]] ^= 1;
				}
			}

			//add to basis if possible
			for (auto it = remainingPositions.begin(); it != remainingPositions.end(); ++it) {
				if (eliminatedVector[*it] != 1) {
					continue;
				}
				//found new basis vector
				exponentBaseMatrix[*it] = exponentVector;
				listingT[*it] = T;
				for (uint64_t j = 0; j < primeBaseSize; ++j) {
					if (j != *it && eliminatedVector[j] == 1) {
						eliminatedVector[j] = 0;
						operationsListVector[*it].push_back(j);
					}
				}
				fillOrder[lastBasisPos] = *it;
				remainingPositions.erase(it);
				++lastBasisPos;
				return 0;
			}

			//calculate factorization
			mpz_class left = T;
			//cout << "factors "<<T<<" ";
			for (uint64_t i = 0; i < primeBaseSize; ++i) {
				if (eliminatedVector[i] != 1) {
					continue;
				}
				//cout << listingT[i]<<" ";
				left = (left * listingT[i]) % n;
				for (uint64_t k = 0; k < primeBaseSize; ++k) {
					exponentVector[k] += exponentBaseMatrix[i][k];
				}
			}
			//cout << endl;
			mpz_class right = 1;
			mpz_class multiplier;
			for (uint64_t k = 0; k < primeBaseSize; ++k) {
				mpz_pow_ui(multiplier.get_mpz_t(),
					mpz_class(primes[k]).get_mpz_t(),
					mpz_class(exponentVector[k] / 2).get_ui());
				right = (right*multiplier) % n;
			}
			if (left != right) {
				//multiplier = gcd(abs(left - right), n);
				left -= right;
				mpz_gcd(multiplier.get_mpz_t(), left.get_mpz_t(), n.get_mpz_t());
				if (abs(multiplier) == 1 || n % multiplier != 0)
					return 0;
				std::cout << "multiplier " << multiplier << endl;
				std::cout << endl << endl;
				return multiplier;
			}
			return 0;
		}
	};

	void step7(mpz_class n, mpz_class& A, mpz_class& T, uint64_t& TSqr, vector<uint64_t>& exponentMatrixRow) {
		if (n % 8 != 1) {
			//if (T % 2 == 1) {
			if (mpz_odd_p(T.get_mpz_t())) {
				exponentMatrixRow[0] = 1;
				TSqr /= 2;
			}
		}
		else {
			mpz_class lowerBound = sqrt(n) + 1;
			mpz_class upperBound = sqrt(n) + A;

			mpz_class residue = n;
			mpz_class modBase = 2;
			mpz_class currentMod = 2;

			mpz_class b = 3;

			mpz_class root1 = 1;
			mpz_class root2 = 3;
			mpz_class root3 = 5;
			mpz_class root4 = 7;
			mpz_class t1 = 1;
			mpz_class t2 = 3;
			mpz_class t3 = 5;
			mpz_class t4 = 7;
			while (!chekNuberInRange(lowerBound, upperBound, root1, currentMod, A)
				&& !chekNuberInRange(lowerBound, upperBound, root2, currentMod, A)
				&& !chekNuberInRange(lowerBound, upperBound, root3, currentMod, A)
				&& !chekNuberInRange(lowerBound, upperBound, root4, currentMod, A)) {
				t1 = root1;
				t2 = root2;
				t3 = root3;
				t4 = root4;
				root1 = liftRootEvenModulo(root1, n, currentMod);
				root2 = liftRootEvenModulo(root2, n, currentMod);
				root3 = liftRootEvenModulo(root3, n, currentMod);
				root4 = liftRootEvenModulo(root4, n, currentMod);
				currentMod *= modBase;
				++b;
			}
			do {
				t1 = root1;
				t2 = root2;
				t3 = root3;
				t4 = root4;
				root1 = liftRootEvenModulo(root1, n, currentMod);
				root2 = liftRootEvenModulo(root2, n, currentMod);
				root3 = liftRootEvenModulo(root3, n, currentMod);
				root4 = liftRootEvenModulo(root4, n, currentMod);
				currentMod *= modBase;
				++b;
			} while (chekNuberInRange(lowerBound, upperBound, root1, currentMod, A)
				|| chekNuberInRange(lowerBound, upperBound, root2, currentMod, A)
				|| chekNuberInRange(lowerBound, upperBound, root3, currentMod, A)
				|| chekNuberInRange(lowerBound, upperBound, root4, currentMod, A));

			--b;
			mpz_class mod;
			mpz_class differ;

			for (uint64_t k = 0; k < b; ++k) {
				differ = abs(T - t1);
				mpz_pow_ui(mod.get_mpz_t(), mpz_class(2).get_mpz_t(), k + 1);
				if (differ % mod == 0) {
					if (exponentMatrixRow[0] < k + 1) {
						exponentMatrixRow[0] = k + 1;
						//step 6
						TSqr /= 2;
					}
				}
				differ = abs(T - t2);
				if (differ % mod == 0) {
					if (exponentMatrixRow[0] < k + 1) {
						exponentMatrixRow[0] = k + 1;
						//step 6
						TSqr /= 2;
					}
				}
				differ = abs(T - t3);
				if (differ % mod == 0) {
					if (exponentMatrixRow[0] < k + 1) {
						exponentMatrixRow[0] = k + 1;
						//step 6
						TSqr /= 2;
					}
				}
				differ = abs(T - t4);
				if (differ % mod == 0) {
					if (exponentMatrixRow[0] < k + 1) {
						exponentMatrixRow[0] = k + 1;
						//step 6
						TSqr /= 2;
					}
				}
			}

		}
	}

	void step456(mpz_class& T, uint64_t& TSqr, vector<uint64_t>& exponentVector, vector<uint64_t>& PrimeBPowers,
		vector<std::pair<mpz_class, mpz_class>>& SqrtmoduloSolutions, vector<uint64_t>& factorBase, uint64_t& factorBaseSize) {
		mpz_class mod;
		mpz_class differ;
		uint64_t b;
		mpz_class t1;
		mpz_class t2;
		uint64_t p;
		for (uint64_t i = 1; i < factorBaseSize; ++i) {
			b = PrimeBPowers[i];
			t1 = SqrtmoduloSolutions[i].first;
			t2 = SqrtmoduloSolutions[i].second;
			p = factorBase[i];

			for (uint64_t k = 0; k < b; ++k) {
				differ = abs(T - t1);
				mpz_pow_ui(mod.get_mpz_t(), mpz_class(p).get_mpz_t(), k + 1);
				//if (differ % (mpz_class)pow(p, k + 1) == 0) {
				if (differ % mod == 0) {
					if (exponentVector[i] < k + 1) {
						exponentVector[i] = k + 1;
						//step 6
						TSqr /= p;
					}
				}
				differ = abs(T - t2);
				//if (differ % (mpz_class)pow(p, k + 1) == 0) {
				if (differ % mod == 0) {
					if (exponentVector[i] < k + 1) {
						exponentVector[i] = k + 1;
						//step 6
						TSqr /= p;
					}
				}
			}
		}
	}

	mpz_class step45678(mpz_class n, mpz_class A, vector<uint64_t>& factorBase) {
		GaussianEliminator Gauss(factorBase);
		cout << "step 4\n";
		uint64_t factorBaseSize = factorBase.size();

		//step 3
		vector<uint64_t> PrimeBPowers(factorBaseSize);
		vector<std::pair<mpz_class, mpz_class>> SqrtmoduloSolutions(factorBaseSize);
		for (uint64_t i = 1; i < factorBaseSize; ++i) {//for each odd prime
			uint64_t p = factorBase[i];
			mpz_class t1 = 0, t2 = 0;
			PrimeBPowers[i] = Step4(n, p, A, t1, t2);
			SqrtmoduloSolutions[i] = { t1,t2 };
		}
		//step 2
		cout << "step 2\n";
		vector<uint64_t> exponentVector(factorBaseSize, 0);
		mpz_class T;
		uint64_t TSqr;
		mpz_class sqrtN = sqrt(n);
		uint64_t printCount = A.get_ui();
		int number = 0;
		float step = 1;
		while (printCount > 10000) {
			printCount /= 10;
			++number;
			step /= 10;
		}
		float curprint = step;
		cout << "Pc " << printCount << endl;

		char buff[100];


		for (uint64_t j = 0; j < A; ++j) {
			if (j% printCount == 0) {
				//printf("\r%.*f%%", number, curprint);
				snprintf(buff, sizeof(buff), "progress \r%.*f%%", number, curprint);
				std::string buffAsStdStr = buff;
				sendProgress(buffAsStdStr);
				curprint += step;
			}
			std::fill(exponentVector.begin(), exponentVector.end(), 0);
			T = sqrtN + j + 1;
			TSqr = mpz_class(T * T - n).get_ui();
			//step 4
			step456(T, TSqr, exponentVector, PrimeBPowers,
				SqrtmoduloSolutions, factorBase, factorBaseSize);
			//step 7
			step7(n, A, T, TSqr, exponentVector);
			//step 8
			if (TSqr != 1)
				continue;
			mpz_class factor = Gauss.addVector(n, T, exponentVector);
			if (factor > 0) {
				//cout << printCount << " " << j << endl;
				return factor;
			}
		}
		return 0;
	}

	void getBound(mpz_class& n, uint64_t& P, uint64_t& A) {
		string numberString = n.get_str();
		mpfr_rnd_t rnd = mpfr_get_default_rounding_mode();
		mpfr_t exp;
		mpfr_init_set_str(exp, numberString.c_str(), 10, rnd);
		mpfr_t logN;
		mpfr_t loglogN;
		mpfr_init(logN);
		mpfr_init(loglogN);
		mpfr_log(logN, exp, rnd);
		mpfr_log(loglogN, logN, rnd);
		mpfr_mul(exp, logN, loglogN, rnd);
		mpfr_sqrt(exp, exp, rnd);
		mpfr_exp(exp, exp, rnd);
		mpfr_rint_ceil(exp, exp, rnd);
		mp_exp_t exponent;
		char* str = mpfr_get_str(NULL, &exponent, 10, 0, exp, rnd);
		numberString = string(str);
		mpfr_free_str(str);
		mpfr_clear(exp);
		P = std::stoull(numberString.substr(0, exponent)) / 10;
		A = P * 10;
	}

	mpz_class QS(mpz_class n, Atkin& atkin) {
		//n is odd
		//step 1
		uint64_t P;
		uint64_t A;
		getBound(n, P, A);
		cout << "A = " << A << endl << "P = " << P << endl;

		if (atkin.primes.back() < P)
			atkin = Atkin(P);
		//step 3
		cout << "step 3\n";
		vector<uint64_t> factorBase;
		uint64_t size = atkin.primes.size();
		factorBase.push_back(2);
		for (uint64_t i = 1; i < size; ++i) {//for each odd prime
			if (atkin.primes[i] > P)
				break;
			//if (LegendreSymbol(n, atkin.primes[i]) == 1)
			if (mpz_legendre(n.get_mpz_t(), mpz_class(atkin.primes[i]).get_mpz_t()) == 1)
				factorBase.push_back(atkin.primes[i]);
		}
		//step 4 5 6 7 8 0
		return step45678(n, A, factorBase);
	}

	int primecheck(mpz_class n) {//0-composite, 1-prime
		if (n == 1 || n == 2)
			return 1;
		if (!mil_rab(n))
			return 0;
		else
			if (!Miller(n))
			{
				return 0;
			}
			else
				return 1;
	}

	uint64_t dumbcheck(mpz_class num, Atkin& atkin) {
		uint64_t size = atkin.primes.size();
		for (uint64_t i = 0; i < size; ++i)
			if (num % atkin.primes[i] == 0)
				return atkin.primes[i];
		return 0;
	}

	vector<string> fatorize(string n)
	{
		vector<string> result;
		Atkin atkin(1000000);
		//string n = "584100001812560001217223";
		mpz_class num;
		mpz_set_str(num.get_mpz_t(), n.c_str(), 10);
		vector<mpz_class> factors;
		vector<mpz_class> numbers = { num };
		while (numbers.size() > 0) {
			num = numbers.back();
			//cout << "factorizing " << num << endl;
			numbers.pop_back();
			if (primecheck(num)) {
				factors.push_back(num);
				continue;
			}
			mpz_class factor;
			factor = dumbcheck(num, atkin);
			if (factor > 0) {
				factors.push_back(factor);
				numbers.push_back(num / factor);
				continue;
			}
			factor = QS(num, atkin);
			if (factor == 0) {
				//cout << "factorization error\n";
				result.push_back("factorization error");
				return result;
			}
			numbers.push_back(factor);
			numbers.push_back(num / factor);

		}
		//cout << "factors\n";
		for (auto it = factors.begin(); it != factors.end(); ++it) {
			mpz_class factor = *it;
			result.push_back(string(mpz_get_str(NULL, 10, factor.get_mpz_t())) + " ");
			//cout << *it << endl;
		}
		//std::cout << endl;
		return result;
	}
}
