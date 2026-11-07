#include "CollieQuickKEYS.hh"

#ifndef SQR
/// Returns the square of its argument
#define SQR(x) ((x) * (x))
#endif

CollieQuickKEYS::CollieQuickKEYS() {
  m_centers = NULL;
  n_fills = 0;
}

CollieQuickKEYS::~CollieQuickKEYS() {
  if (m_centers != NULL)
    delete[] m_centers;
}

// overrides
int CollieQuickKEYS::Book(int id, const char* name, int bins, double left, double right,
                          CollieHistogramContainer* container) {
  CollieHistogram::Book(id, name, bins, left, right, container);
  if (m_centers != NULL)
    delete[] m_centers;
  m_centers = new double[bins];
  for (int i = 0; i < bins; i++)
    m_centers[i] = left + (right - left) / double(bins) * double(i);
  return true;
}

int CollieQuickKEYS::Book(const char* name, int bins, double left, double right,
                          CollieHistogramContainer* container) {
  CollieHistogram::Book(name, bins, left, right, container);
  if (m_centers != NULL)
    delete[] m_centers;
  m_centers = new double[bins];
  for (int i = 0; i < bins; i++)
    m_centers[i] = left + (right - left) / double(bins) * double(i);
  return true;
}

void CollieQuickKEYS::Fill(double x, double w) {
  if (m_centers != NULL && x >= m_start && x < m_end) {
    int bin = int((x - m_start) / (m_end - m_start) * double(n_bins));
    m_centers[bin] = (m_centers[bin] * p_Histogram[bin] + x * w) / (p_Histogram[bin] + w);
    n_fills++;
  }
  CollieHistogram::Fill(x, w);
}

void CollieQuickKEYS::Fill(int i, double w) {
  if (m_centers != NULL && i >= 0 && i < n_bins) {
    double center = (double(i) + 0.5) * (m_end - m_start) / double(n_bins) + m_start;
    m_centers[i] = (m_centers[i] * p_Histogram[i] + center * w) / (p_Histogram[i] + w);
    n_fills++;
  }
  CollieHistogram::Fill(i, w);
}

void CollieQuickKEYS::Clear() {
  for (int i = 0; i < n_bins; i++)
    m_centers[i] = 0;
  n_fills = 0;
  CollieHistogram::Clear();
}

void CollieQuickKEYS::Copy(double* copy, int bins_per_bin) {
  int finalbins = (n_bins + bins_per_bin - 1) / bins_per_bin;
  for (int i = 0; i < finalbins; i++) {
    copy[i] = p_Histogram[i * bins_per_bin];
    for (int j = 1; j < bins_per_bin; j++)
      copy[i] += p_Histogram[i * bins_per_bin + j];
  }
}

#define T(j) (((double)(j) + 0.5f) * binw + m_start)
#define T1(j) (key->x[j])
#define T2(j) ((j < 0 || j >= n_bins) ? T(j) : T1(j))
#define Gaus(x) ((fabs(x) > 10.0) ? 0.0 : exp(-((x) * (x)) / 2.0))

void CollieQuickKEYS::Smooth(double m_alpha, bool refl_f0) {

  int i, j;
  double *pass1, *pass2;
  double h_i, h, sw, xave, sigma;
  double binw;

  /* smoothing is now good */
  binw = (m_end - m_start) / ((double)n_bins);

  if (n_fills == 0) {
    for (i = 0; i < n_bins; i++) {
      p_Histogram[i] = 0;
      p_Fills[i] = 1;
    }
    return;
  }

  sw = 0;
  xave = 0;
  sigma = 0;

  /* pass 1 creates the "f" distribution */
  pass1 = new double[n_bins];
  for (i = 0; i < n_bins; i++) {
    pass1[i] = 0;
    sw += p_Histogram[i];
    xave += m_centers[i] * p_Histogram[i];
  }
  xave /= sw;

  if (n_fills == 1) {
    for (i = 0; i < n_bins; i++) {
      p_Histogram[i] = sw / (double)n_bins;
      p_Fills[i] = 1;
    }
    return;
  }

  for (i = 0; i < n_bins; i++)
    sigma += SQR(m_centers[i] - xave) * p_Histogram[i];
  sigma = sqrt(sigma / sw);
  h = 1.059223841 * pow(n_fills, -0.2) * sigma;
  if (h < binw)
    h = binw;

  /*  add up the gaussians */
  for (i = 0; i < n_bins; i++)
    for (j = 0; j < n_bins; j++)
      if (p_Histogram[i] > 0) {
        pass1[j] += p_Histogram[i] * Gaus((m_centers[j] - m_centers[i]) / h);
      }

  /* do the reflections, if requested */
  if (refl_f0) {
    for (i = 0; i < n_bins; i++) {
      for (j = 0; j < n_bins; j++)
        if (p_Histogram[j] > 0) {
          /* lower reflection is at -i-1 */
          pass1[i] += p_Histogram[j] * Gaus((m_centers[j] - T(-i - 1)) / h);
          /* upper reflection is at 2*n_bins-i-1 (2*n_bins-n_bins+1-1=n_bins) */
          pass1[i] += p_Histogram[j] * Gaus((m_centers[j] - T(2 * n_bins - i - 1)) / h);
        }
    }
  }

  /* normalize pass1 */
  for (i = 0; i < n_bins; i++) {
    pass1[i] /= sqrt(2.0 * M_PI) * sw * h;
  }
  h = 1.059223841 * pow(n_fills, -0.2) * sqrt(sigma) / (2.0 * sqrt(3.0)) * m_alpha;

  /* pass 2 creates the final histogram */
  pass2 = new double[n_bins * 3];
  for (i = -n_bins; i < 2 * n_bins; i++) {
    pass2[i + n_bins] = 0;
    for (j = 0; j < n_bins; j++)
      if (p_Histogram[j] > 0 && pass1[j] > 0) {
        h_i = h / sqrt(pass1[j]);
        if (h_i < h * sigma / 10.0)
          h_i = h * sigma / 10.0;
        if (h_i < binw)
          h_i = binw;
        pass2[i + n_bins] += p_Histogram[j] * Gaus((T(i) - m_centers[j]) / h_i) / h_i;
      }
  }

  h = 0; /* double use of h -- ignore name */
  /* normalize pass2 */
  for (i = 0; i < n_bins; i++) {
    pass2[i + n_bins] += pass2[n_bins - i - 1] + pass2[n_bins * 3 - i - 1];
    pass2[i + n_bins] /= 2.0 * M_PI;
    h += pass2[i + n_bins];
  }

  /* normal into the final histogram */
  for (i = 0; i < n_bins; i++) {
    p_Histogram[i] = pass2[i + n_bins] * sw / h;
    p_Fills[i] = 1;
  }

  delete[] pass1;
  delete[] pass2;
}
