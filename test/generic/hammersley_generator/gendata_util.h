int get_prime(int n);
void init_fib_rnd(double *rnd_seed, int *fib_i, int *fib_j, int *i_ctrl);
double fib_rnd(double *x, double *rnd_seed, int *fib_i, int *fib_j,
               int *i_ctrl);
void set_random_charge(int n, double *q);
void set_negative_charge(int n, double *q);
void set_positive_charge(int n, double *q);
void shift_scale(int d, int n, double *x, double shift, double scale);
void auto_shift_scale(int d, int n, double *x);
