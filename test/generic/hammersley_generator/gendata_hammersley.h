void gendata_hammersley_cube(int n, double **x, int *n_total);
void gendata_hammersley_ball(int n, double **x, int *n_total);
void gendata_hammersley_ball_neg_charge(int n, double **x, int *n_total);
void gendata_hammersley_sphere(int n, double **x, int *n_total);
void gendata_hammersley_circle(int n, double **x, int *n_total);
void gendata_hammersley_square(int n, double **x, int *n_total);
void gendata_halton_ellipsoid(int n, double **x, double a, double b,
                              double c, int *n_total);
void gendata_halton_cylinder(int n, double **x, double r_div_len, int *n_total);
void gendata_hammersley_two_balls(int n, double **x, double n_2_rel_n_1,
                                  double distance_rel_r_1, int *n_total);
