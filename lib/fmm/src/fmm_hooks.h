
#ifndef __FMM_HOOKS_H__
#define __FMM_HOOKS_H__


typedef struct _fcs_fmm_hooks_t
{
  int demo_value;
  double t_near, t_near_d, t_far, t_far_d;

} fcs_fmm_hooks_t;

#define FCS_FMM_HOOKS_NULL  NULL


fcs_fmm_hooks_t *fcs_fmm_hooks_create();
void fcs_fmm_hooks_destroy(fcs_fmm_hooks_t *hooks);

void fcs_fmm_hooks_near_start(fcs_fmm_hooks_t *hooks);
void fcs_fmm_hooks_near_stop(fcs_fmm_hooks_t *hooks);
void fcs_fmm_hooks_far_start(fcs_fmm_hooks_t *hooks);
void fcs_fmm_hooks_far_stop(fcs_fmm_hooks_t *hooks);

void fcs_fmm_hooks_print(fcs_fmm_hooks_t *hooks);


#endif /* __FMM_HOOKS_H__ */
