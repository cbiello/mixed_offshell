#include "handyG_cpp.h"

extern "C" void g_weight1_fortran (const double * re1, const double * im1, const double * reX, const double * imX, double * re_ans, double* im_ans);
extern "C" void g_weight2_fortran (const double * re1, const double * im1, const double * re2, const double * im2, const double * reX, const double * imX, double * re_ans, double* im_ans);
extern "C" void g_weight3_fortran (const double * re1, const double * im1, const double * re2, const double * im2, const double * re3, const double * im3, const double * reX, const double * imX, double * re_ans, double* im_ans);
extern "C" void g_weight4_fortran (const double * re1, const double * im1, const double * re2, const double * im2, const double * re3, const double * im3, const double * re4, const double * im4, const double * reX, const double * imX, double * re_ans, double* im_ans);


namespace handyG
{
	std::complex<double> G (const std::complex<double> a1, const std::complex<double> x)
	{
		const double a1_re = a1.real();
		const double a1_im = a1.imag();
		const double x_re = x.real();
		const double x_im = x.imag();
		double ans_re;
		double ans_im;
		g_weight1_fortran(&a1_re, &a1_im, &x_re, &x_im, &ans_re, &ans_im);
		return std::complex<double>(ans_re, ans_im);
	}
	std::complex<double> G (const std::complex<double> a1, const std::complex<double> a2, const std::complex<double> x)
	{
		const double a1_re = a1.real();
		const double a1_im = a1.imag();
		const double a2_re = a2.real();
		const double a2_im = a2.imag();
		const double x_re = x.real();
		const double x_im = x.imag();
		double ans_re;
		double ans_im;
		g_weight2_fortran(&a1_re, &a1_im, &a2_re, &a2_im, &x_re, &x_im, &ans_re, &ans_im);
		return std::complex<double>(ans_re, ans_im);
	}
	std::complex<double> G (const std::complex<double> a1, const std::complex<double> a2, const std::complex<double> a3, const std::complex<double> x)
	{
		const double a1_re = a1.real();
		const double a1_im = a1.imag();
		const double a2_re = a2.real();
		const double a2_im = a2.imag();
		const double a3_re = a3.real();
		const double a3_im = a3.imag();
		const double x_re = x.real();
		const double x_im = x.imag();
		double ans_re;
		double ans_im;
		g_weight3_fortran(&a1_re, &a1_im, &a2_re, &a2_im, &a3_re, &a3_im, &x_re, &x_im, &ans_re, &ans_im);
		return std::complex<double>(ans_re, ans_im);
	}
	std::complex<double> G (const std::complex<double> a1, const std::complex<double> a2, const std::complex<double> a3, const std::complex<double> a4, const std::complex<double> x)
	{
		const double a1_re = a1.real();
		const double a1_im = a1.imag();
		const double a2_re = a2.real();
		const double a2_im = a2.imag();
		const double a3_re = a3.real();
		const double a3_im = a3.imag();
		const double a4_re = a4.real();
		const double a4_im = a4.imag();
		const double x_re = x.real();
		const double x_im = x.imag();
		double ans_re;
		double ans_im;
		g_weight4_fortran(&a1_re, &a1_im, &a2_re, &a2_im, &a3_re, &a3_im, &a4_re, &a4_im, &x_re, &x_im, &ans_re, &ans_im);
		return std::complex<double>(ans_re, ans_im);
	}
}
