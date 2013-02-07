void* ptroffset_lli(long long *ptr, int offset)
{
return ptr + offset;
}
void* ptroffset_li(long *ptr, int offset)
{
return ptr + offset;
}
void* ptroffset_pi(void *ptr, int offset)
{
return (void *) ((char *)ptr + offset);
}
void* ptroffset_ii(int *ptr, int offset)
{
return ptr + offset;
}
void* ptroffset_fi(float *ptr, int offset)
{
return ptr + offset;
}
void* ptroffset_di(double *ptr, int offset)
{
return ptr + offset;
}

void* ptroffset_llll(long long *ptr, long long offset)
{
return ptr + offset;
}
void* ptroffset_lll(long *ptr, long long offset)
{
return ptr + offset;
}
void* ptroffset_pll(void *ptr, long long offset)
{
return (void *) ((char *)ptr + offset);
}
void* ptroffset_ill(int *ptr, long long offset)
{
return ptr + offset;
}
void* ptroffset_fll(float *ptr, long long offset)
{
return ptr + offset;
}
void* ptroffset_dll(double *ptr, long long offset)
{
return ptr + offset;
}
