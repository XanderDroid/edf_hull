#include "edf_hull.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IND (cur_points->vec_sel[i])

void edf_set_zero(edf_points_t *cur_points) {
  cur_points->alloc_points = 0;
  cur_points->t0 = NULL;
  cur_points->t1 = NULL;
  cur_points->vec_p = NULL;
  cur_points->num_sel = 0;
  cur_points->vec_sel = NULL;
}

void edf_create_points(const ts_t *cur_task_set, edf_points_t *cur_points) {
#define N (cur_task_set->num)
#define T(i) (cur_task_set->per[i])
#define D(i) (cur_task_set->dl[i])
#define O(i) (cur_task_set->phi[i])

  unsigned long *max_job;
  unsigned long i, j, k, s_i, s_j;
  int reserve_p;
  int num_deadline;
  int point_count;
  double cur_t, cur_s, coef;
  double big_enough;

  /* store the number of tasks */
  cur_points->num_tasks = N;

  /* Reserve one point for sum_i U_i <= 1, N points for the Ci >= 0 */
  reserve_p = N + 1;

  /* Limit until which computing the points + number of points */
  if (cur_task_set->has_phi) {
    big_enough =
        2 * cur_task_set->h_per + cur_task_set->max_d + cur_task_set->max_o;
    /* the number of points (overestimating it) */
    cur_points->num_points = reserve_p + num_deadline * num_deadline;
  } else {
    big_enough = cur_task_set->h_per + cur_task_set->max_d;
    /* the number of points */
    cur_points->num_points = reserve_p + num_deadline;
  }

  /*
   * Maximum job of the tasks. Also compute the sum of max_job[i],
   * because this is also the number of considered deadlines.
   */
  max_job = malloc(sizeof(*max_job) * N);
  num_deadline = 0;
  for (i = 0; i < N; i++) {
    max_job[i] = (unsigned long)floor((big_enough - O(i) - D(i)) / T(i));
    num_deadline += max_job[i] + 1;
  }

  /* allocate the t0, t1 data structure for the points */
  if (cur_points->num_points > cur_points->alloc_points) {
    cur_points->t0 =
      (double *)realloc(cur_points->t0, sizeof(double) * cur_points->num_points);
    cur_points->t1 =
      (double *)realloc(cur_points->t1, sizeof(double) * cur_points->num_points);
    cur_points->vec_p =
      (double *)realloc(cur_points->vec_p, sizeof(double) * cur_points->num_points * N);
    cur_points->alloc_points = cur_points->num_points;
  }

  /* Start computing the points. */
  /* initialize the point counter */
  point_count = 0;

  /* store the -Ci <= 0 constraint */
  for (i = 0; i < N; i++) {
    cur_points->t0[point_count] = 0;
    cur_points->t1[point_count] = 0;
    for (j = 0; j < N; j++) {
      cur_points->vec_p[point_count * N + j] = ((i == j) ? -1 : 0);
    }
    
    /* increment the point counter */
    point_count++;
  }
  
  /* store the sum Ui <= 1 constraint */
  cur_points->t0[point_count] = 0;
  cur_points->t1[point_count] = cur_task_set->h_per;
  for (i = 0; i < N; i++) {
    cur_points->vec_p[point_count * N + i] = (int)(cur_task_set->h_per / T(i));
  }
  point_count++;

  /* if we have offset the points can be really many */
  if (cur_task_set->has_phi) {
    /* compute the points corresponding to the (start time, deadline) pairs */
    /* loop on the absolute deadlines */
    for (i = 0; i < N; i++) {
      for (j = 0; j <= max_job[i]; j++) {
        /* t1 is equal to the j-th deadline of task i */
        cur_t = O(i) + j * T(i) + D(i);

        /* loop on the start times */
        for (s_i = 0; s_i < N; s_i++) {
          for (s_j = 0; s_j <= max_job[s_i]; s_j++) {
            /* the current start time */
            cur_s = O(s_i) + s_j * T(s_i);
            if (cur_s >= cur_t) /* start time is ahead the deadline */
              break;
            cur_points->t1[point_count] = cur_t;
            cur_points->t0[point_count] = cur_s;
            for (k = 0; k < N; k++) {
              /* there may be numerical problems... */
              coef = floor((cur_t - O(k) - D(k)) / T(k)) -
                     ceil((cur_s - O(k)) / T(k)) + 1;
              /* must be non-negative */
              if (coef < 0)
                coef = 0;
              /* store it */
              cur_points->vec_p[point_count * N + k] = coef;
            }
            /* increment the point counter */
            point_count++;
          }
        }
      }
    }
    /* update the num_points */
    cur_points->num_points = point_count;
  } else { /* no offset here */

    /* the number of points */
    cur_points->num_points = reserve_p + num_deadline;

    if (cur_points->num_points > cur_points->alloc_points) {
      /* allocate the t0, t1 data structure for the points */
      cur_points->t0 = (double *)realloc(
          cur_points->t0, sizeof(double) * cur_points->num_points);
      cur_points->t1 = (double *)realloc(
          cur_points->t1, sizeof(double) * cur_points->num_points);
      cur_points->vec_p = (double *)realloc(
          cur_points->vec_p, sizeof(double) * cur_points->num_points * N);
      cur_points->alloc_points = cur_points->num_points;
    }

    /* initialize the point counter */
    point_count = 0;

    /* store the -Ci <= 0 constraint */
    for (i = 0; i < N; i++) {
      cur_points->t0[point_count] = 0;
      cur_points->t1[point_count] = 0;
      for (j = 0; j < N; j++) {
        cur_points->vec_p[point_count * N + j] = ((i == j) ? -1 : 0);
      }

      /* increment the point counter */
      point_count++;
    }

    /* store the sum Ui <= 1 constraint */
    cur_points->t0[point_count] = 0;
    cur_points->t1[point_count] = cur_task_set->h_per;
    for (i = 0; i < N; i++) {
      cur_points->vec_p[point_count * N + i] = (int)(cur_task_set->h_per / T(i));
    }
    point_count++;

    /* compute the points corresponding to deadlines */
    for (i = 0; i < N; i++) {
      for (j = 0; j <= max_job[i]; j++) {
        /* t0 is zero, when no offset */
        cur_points->t0[point_count] = 0;
        /* t1 is equal to the j-th deadline of task i */
        cur_points->t1[point_count] = cur_t = j * T(i) + D(i);
        for (k = 0; k < N; k++) {
          /* there may be numerical problems... */
          coef = floor((cur_t - D(k)) / T(k)) + 1;
          /* must be non-negative */
          if (coef < 0)
            coef = 0;
          /* store it */
          cur_points->vec_p[point_count * N + k] = coef;
        }

        /* increment the point counter */
        point_count++;
      }
    }
  }

  /* free the max_job which is only a local data structure */
  free(max_job);
#undef N
#undef T
#undef D
#undef O
}

void edf_print_constraints_C(const edf_points_t *cur_points) {
  int i, j;
  printf("\nMinimal set of constraints are written in the following form\n\n");
  printf("    eta_1*C_1 +eta_2*C_2 + ... +eta_N*C_N <= t_1-t_0\n\n");
  for (j = 0; j < cur_points->num_tasks; j++) {
    printf("eta_%d\t", j + 1);
  }
  printf("t1\tt0\n");
  for (i = 0; i < cur_points->num_sel; i++) {
    for (j = 0; j < cur_points->num_tasks; j++) {
      printf("%.4f\t", cur_points->vec_p[IND * cur_points->num_tasks + j]);
    }
    printf("%.0f\t%.0f\n", cur_points->t1[IND], cur_points->t0[IND]);
  }
}

double edf_linprog_points(edf_points_t *cur_points) {
  unsigned int i, j;
  double total_t;

  printf("@Alessandro: FIXME START\n");
  printf("Istruzioni:\n");
  printf("1. nella funzione edf_linprog_points trovi questo codice\n");
  printf("\n2. tutti i vincoli sono scritti nel seguente formato\n");
  printf("\ta_1*U_1 +a_2*U_2 + ... +a_N*U_N <= RHS\n");
  printf("\n3. vincoli -Ci <= 0 (vedi te se serve aggiungerli esplicitamente)\n");
  for (i = 0; i < cur_points->num_tasks; i++) {
    for (j = 0; j < cur_points->num_tasks; j++) {
      printf("%.4f\t", cur_points->vec_p[i * cur_points->num_tasks + j]);
    }
    printf("\t\t%6f\n", cur_points->t1[i]-cur_points->t0[i]);   /* RHS */
  }
  printf("\n4. vincolo sum Ci <= 1 (primo vincolo sicuro)\n");
  for (j = 0; j < cur_points->num_tasks; j++) {
    printf("%.4f\t", cur_points->vec_p[i * cur_points->num_tasks + j]);
  }
  printf("\t\t%6f\n", cur_points->t1[i]-cur_points->t0[i]);   /* RHS */
  i++;    /* Move to next constraint */
  printf("\n5. seguono i vincoli delle deadline che vanno aggiungi uno alla volta\n");
  for (/* i */; i < cur_points->num_points; i++) {
    for (j = 0; j < cur_points->num_tasks; j++) {
      printf("%.4f\t", cur_points->vec_p[i * cur_points->num_tasks + j]);
    }
    printf("\t\t%6f\n", cur_points->t1[i]-cur_points->t0[i]);   /* RHS */
  }
  printf("\n6. alla fine, bisogna aggiornare:\n");
  printf("\tcur_points->num_sel con il numero di vincoli selezionati\n");
  printf("\tcur_points->vec_sel array con indici punti selezionati\n");
  /* FIXME: place-holder, just selecting all points */
  cur_points->num_sel = cur_points->num_points;
  cur_points->vec_sel = malloc(sizeof(*cur_points->vec_sel)
			       *cur_points->num_sel);
  for(i=0; i<cur_points->num_sel; i++) {
    cur_points->vec_sel[i]=i;
  }
  printf("@Alessandro: FIXME END\n");
  
  
  total_t = -1; /* FIXME: should be amount of consumed time  */

  return total_t;
}

void edf_free_points(edf_points_t *cur_points) {
  free(cur_points->t0);
  free(cur_points->t1);
  free(cur_points->vec_p);
  free(cur_points->vec_sel);
}

void edf_print_constraints_U(const ts_t *cur_task_set,
                             const edf_points_t *cur_points) {
  int i, j;
  printf("\nOr, alternatively, minimal set of constraints can also be written "
         "as\n\n");
  printf("    a_1*U_1 +a_2*U_2 + ... +a_N*U_N <= t_1-t_0\n\n");
  for (j = 0; j < cur_points->num_tasks; j++) {
    printf("a_%d\t", j + 1);
  }
  printf("t1\tt0\n");

#define T(j) (cur_task_set->per[j])
  for (i = 0; i < cur_points->num_sel; i++) {
    for (j = 0; j < cur_points->num_tasks; j++) {
      if (i < cur_points->num_tasks) {
	/* Special case of -Ci <= 0 constraint */
        printf("%.4f\t", cur_points->vec_p[IND * cur_points->num_tasks + j]);
      } else {
	/* Other constraints: the 1st one is tot util<=1 */
        printf("%.4f\t",
               (cur_points->vec_p[IND * cur_points->num_tasks + j] * T(j)) /
                   cur_points->t1[IND]);
      }
    }
    if (i < cur_points->num_tasks) {
      printf("0\t%.0f\n", cur_points->t0[IND]);
    } else {
      printf("1\t%.0f\n", cur_points->t0[IND]);
    }

#undef T
  }
}

void edf_print_points(const edf_points_t *cur_points) {
  int i, j;

  printf("Number of tasks:\t%d\n", cur_points->num_tasks);
  printf("Number of points:\t%d\n", cur_points->num_points);
  printf("Points:\n");
  for (i = 0; i < cur_points->num_points; i++) {
    printf("t0=%.3f\tt1=%.3f", cur_points->t0[i], cur_points->t1[i]);
    for (j = 0; j < cur_points->num_tasks; j++) {
      printf("\t%.3f", cur_points->vec_p[i * cur_points->num_tasks + j]);
    }
    printf("\n");
  }
}
