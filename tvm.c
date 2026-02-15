#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double n, i, PV, PMT, FV; 

double fn(double n) {
    double t = pow(1.0 + i, n);
    return PV * t + PMT * (t - 1.0) / i + FV;
}

double fi(double i) {
    double t = pow(1.0 + i, n);
    return PV * t + PMT * (t - 1.0) / i + FV;
}

double fn_prime(double n) {
   return log(1.0 + i) * (PV + PMT / i) * pow(1.0 + i, n);
}

double fi_prime(double i) {
    double t = pow(1.0 + i, n);
    double t_minus1 = pow(1.0 + i, n - 1.0);

    return PV * n * t_minus1
         + PMT * ((i * n * t_minus1 - (t - 1.0)) / (i * i));
}

double newton_raphson_n(int line_number) {
    double x= 360.0; 
    double delta = 0.0;
    
    for (int k = 0; k <= 10000; k++) {
        delta = -fn(x) / fn_prime(x);
        x += delta;

        if (fabs(delta) < 1e-8) {
            return ceil(x);
        }
    }

    fprintf(stderr, "line %d: solver did not converge\n", line_number);
    exit(1);
}

double newton_raphson_i(int line_number) {
    double guesses[] = { 0.0025, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0 };
    int num_guesses = (int)(sizeof(guesses) / sizeof(guesses[0]));

    for (int g = 0; g < num_guesses; g++) {
        double x = guesses[g];

    for (int k = 0; k <= 10000; k++) {
        double fx = fi(x);
        double fpx = fi_prime(x);

        if (fx != fx || fpx != fpx || fpx == 0.0) return -NAN;

        double delta = -fx / fpx;
        if (delta != delta) return -NAN;
        
        x += delta;
        if (x <= 0) x = 1e-7;
        if (fabs(delta) < 1e-8) 
            return x;
        }
    }
}
    fprintf(stderr, "line %d: solver did not converge\n", line_number);
    exit(1);

void check_i(int line_number) {
    if (i == 0.0) {
        fprintf(stderr, "line %d: i must not be zero\n", line_number);
        exit(1);
    }
}


void tvm_set_variable(char *name, double value, int line_number) {
    if (strcmp(name, "n") == 0) {
        if (value <= 0 || value != ceil(value)) {
            fprintf(stderr, "line %d: n must be a positive integer\n", line_number);
            exit(1);
    }
    n = value;
    return;
}
    if (strcmp(name, "i") == 0 && value <= 0) {
        fprintf(stderr, "line %d: i must be positive\n", line_number);
        exit(1);
    }
    else if (strcmp(name, "i") == 0) {
        i = value;
    }
    else if (strcmp(name, "PV") == 0) {
        PV = value;
    }
    else if (strcmp(name, "PMT") == 0) {
        PMT = value;
    }
    else if (strcmp(name, "FV") == 0) {
        FV = value;
    }


    else {
        fprintf(stderr, "line %d: unknown variable\n", line_number);
        exit(1);
    }
}


void tvm_compute_variable(char *name,int line_number) {
    if (strcmp(name, "n") == 0) {
        check_i(line_number);
        n = newton_raphson_n(line_number);
        printf("n = %.0f\n", n);
    }
    else if (strcmp(name, "i") == 0) {
        i = newton_raphson_i(line_number);
        printf("i = %.6f\n", i);
}
    else if (strcmp(name, "PV") == 0) {
        check_i(line_number);
        PV = - ( FV + PMT * (pow(1 + i, n) - 1) / i ) / pow(1 + i, n);
        printf("PV = %.2f\n", PV);
    }

    else if (strcmp(name, "PMT") == 0) {
        check_i(line_number);
        PMT = - ( PV * pow(1 + i, n) + FV ) / ( (pow(1 + i, n) - 1) / i );
        printf("PMT = %.2f\n", PMT);
    }
    else if (strcmp(name, "FV") == 0) {
        check_i(line_number);
        FV = -( PV * pow(1 + i, n) + PMT * (pow(1 + i, n) - 1) / i );
        printf("FV = %.2f\n", FV);
    }
    else {
        fprintf(stderr, "line %d: unknown variable\n", line_number);
        exit(1);
    }
}

void tvm_clear(void) {
    
    n = 0.0; 
    i = 0.0; 
    PV = 0.0; 
    PMT = 0.0;
    FV = 0.0; 
}

void tvm_process_command(char *command, int line_number) {
    char varname[10];
    double value;  

    if (sscanf(command, "set %s %lf", varname, &value) == 2 ) {
        tvm_set_variable(varname, value, line_number);
        return; 
    }

    if (sscanf(command, "compute %s", varname) == 1) {
        tvm_compute_variable(varname, line_number);
        return;
    }

    if (strcmp(command, "clear")==0) {
    tvm_clear();
    return; 
    }
}

void truncate_at_newline(char *s) {
    int i = 0;

    while (s[i] != '\0' && s[i] != '\n') ++i;

    s[i] = '\0';
}

int main(void) {
    char command[40];

    int line_number = 0;

    tvm_clear();

    while (1) {
        char *res = fgets(command, sizeof(command), stdin);

        if (res == NULL) break;

        ++line_number;
        truncate_at_newline(command);
        tvm_process_command(command, line_number);
    }

    if (ferror(stdin)) {
        fprintf(stderr, "tvm:  Error reading input\n");
        exit(1);
    }

    return 0;
    }
