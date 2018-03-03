#include <cstdint>
#include <vector>
#include <set>
#include <algorithm>
#include <mpirxx.h>
using namespace std;
class Atkin {
private:
	vector<bool> is_prime;	
public:
	vector<uint64_t> primes = { 2,3,5 };

	Atkin(uint64_t limit) {
		uint64_t sqr_lim;
		is_prime.resize(limit + 1, false);
		uint64_t x2, y2;
		uint64_t i, j;
		uint64_t n;
		// Инициализация решета
			sqr_lim = sqrt(limit);
		for (i = 0; i <= limit; i++) is_prime[i] = false;
		is_prime[2] = true;
		is_prime[3] = true;

		// Предположительно простые — это целые с нечётным числом
		// представлений в данных квадратных формах.
		// x2 и y2 — это квадраты i и j (оптимизация).
		x2 = 0;
		for (i = 1; i <= sqr_lim; i++) {
			x2 += 2 * i - 1;
			y2 = 0;
			for (j = 1; j <= sqr_lim; j++) {
				y2 += 2 * j - 1;

				n = 4 * x2 + y2;
				if ((n <= limit) && (n % 12 == 1 || n % 12 == 5))
					is_prime[n] = !is_prime[n];

				// n = 3 * x2 + y2; 
				n -= x2; // Оптимизация
				if ((n <= limit) && (n % 12 == 7))
					is_prime[n] = !is_prime[n];

				// n = 3 * x2 - y2;
				n -= 2 * y2; // Оптимизация
				if ((i > j) && (n <= limit) && (n % 12 == 11))
					is_prime[n] = !is_prime[n];
			}
		}

		// Отсеиваем кратные квадратам простых чисел в интервале [5, sqrt(limit)].
		// (основной этап не может их отсеять)
		for (i = 5; i <= sqr_lim; i++) {
			if (is_prime[i]) {
				n = i * i;
				for (j = n; j <= limit; j += n) {
					is_prime[j] = false;
				}
			}
		}

		// Вывод списка простых чисел в консоль.
		for (i = 6; i <= limit; i++) {  // добавлена проверка делимости на 3 и 5. В оригинальной версии алгоритма потребности в ней нет.
			if ((is_prime[i]) && (i % 3 != 0) && (i % 5 != 0)) {
				primes.push_back(i);
			}
		}
	}
};
