#if !defined(TIMER_HEADER_)
#define TIMER_HEADER_

/** Start timing. */
void tic (void);

/** Return seconds since last tic. */
double toc (void);

/** Macro to time a block. */
#define TIME(timevar, what) do { tic (); what; timevar = toc(); } while (0)

#endif /* TIMER_HEADER_ */
