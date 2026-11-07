#include "Loader.hh"
#include <string.h>

Loader::Loader() {
  nscales = 0;
  memset(scales, 0, sizeof(scales));
  m_filename = NULL;
  m_options = NULL;
}

Loader::~Loader() {
  if (m_filename != NULL)
    delete[] m_filename;
  if (m_options != NULL)
    delete[] m_options;
  m_filename = NULL;
  m_options = NULL;
}

bool Loader::open(const char* filename, const char* options) {
  if (filename == NULL || filename[0] == '\0')
    return false;
  delete[] m_filename;
  delete[] m_options;
  m_filename = NULL;
  m_options = NULL;
  m_filename = new char[strlen(filename) + 1];
  strcpy(m_filename, filename);
  if (options != NULL) {
    m_options = new char[strlen(options) + 1];
    strcpy(m_options, options);
  }
  return true;
}
