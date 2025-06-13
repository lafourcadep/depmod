#include "parsers.h"
#include <cassert>
#include <filesystem>

namespace io {

/* ------------------------------------------------------------------------- */

TextFileHandler::TextFileHandler(const std::string& filepath) : m_path(filepath) {}

ASCIIFileHandler::ASCIIFileHandler(const std::string& path, FileMode mode) : TextFileHandler(path) {
  m_fptr = std::fopen(m_path.c_str(), convert_mode_to_char(mode));
  if (m_fptr == nullptr) {
    PANIC("Could not open file");
  }
}

ASCIIFileHandler::~ASCIIFileHandler() {
  if (m_fptr != nullptr)
    std::fclose(m_fptr);
}

void ASCIIFileHandler::clear() noexcept { std::clearerr(m_fptr); }

void ASCIIFileHandler::seek(uint64_t cursor) {
  static_assert(sizeof(uint64_t) == sizeof(off64_t));

  auto status = fseeko(m_fptr, static_cast<off64_t>(cursor), SEEK_SET);
  if (status != 0) {
    PANIC("Fail to seek");
  }
}

size_t ASCIIFileHandler::read(char* data, size_t size) {
  size_t result = std::fread(data, 1, size, m_fptr);
  if (std::ferror(m_fptr) != 0) {
    PANIC("Fail to read file");
  }
  return result;
}

#ifdef USE_ZLIB

GzipFileHandler::GzipFileHandler(const std::string& path, FileMode mode) : TextFileHandler(path) {
  const char* openmode;
  switch (mode) {
  case FileMode::READ:
    openmode = "rb";
    break;
  case FileMode::APPEND:
  case FileMode::WRITE:
    PANIC("FileMode::WRITE is not supported");
    break;
  }

  m_fptr = gzopen64(m_path.c_str(), openmode);
  if (m_fptr == nullptr) {
    PANIC("Unable to open file: ", m_path.c_str());
  }
}

GzipFileHandler::~GzipFileHandler() {
  if (m_fptr != nullptr)
    gzclose(m_fptr);
}

void GzipFileHandler::clear() noexcept { gzclearerr(m_fptr); }

void GzipFileHandler::seek(uint64_t cursor) {
  static_assert(sizeof(uint64_t) == sizeof(z_off64_t));

  auto status = gzseek64(m_fptr, static_cast<z_off64_t>(cursor), SEEK_SET);
  if (status == -1) {
    const char* message = gz_error();
    PANIC("Error while decompressing gzip file : %s", message);
  }
}

size_t GzipFileHandler::read(char* buffer, size_t size) {
  int read = gzread(m_fptr, buffer, safe_cast(size));
  const char* error = gz_error();
  if (read == -1 || error != nullptr) {
    PANIC("Error while reading gzip file : %s", error);
  }
  return static_cast<size_t>(read);
}

const char* GzipFileHandler::gz_error(void) const {
  int status = Z_OK;
  const char* message = gzerror(m_fptr, &status);
  return status == Z_OK ? nullptr : message;
}

unsigned GzipFileHandler::safe_cast(size_t value) {
  constexpr size_t max = std::numeric_limits<unsigned>::max();
  if (value > max) {
    PANIC("%d is too big for unsigned in call to zlib function", value);
  }
  return static_cast<unsigned>(value);
}

#endif

#ifdef USE_LZMA

size_t XzFileHandler::safe_cast(size_t size) {
  constexpr size_t max = std::numeric_limits<size_t>::max();
  if (size > max) {
    PANIC("%ld is too big for unsigned in call to zlib function", size);
  }
  return static_cast<size_t>(size);
}

void XzFileHandler::check_lzma_ret(lzma_ret ret) {
  int retnum = static_cast<int>(ret);

  switch (ret) {
  case LZMA_OK:
  case LZMA_STREAM_END:
    return;
  case LZMA_GET_CHECK:
  case LZMA_NO_CHECK:
    PANIC("lzma: no integrity check.", retnum);
  case LZMA_MEM_ERROR:
  case LZMA_MEMLIMIT_ERROR:
    PANIC("lzma: failed to allocated memory (err: %d)", retnum);
  case LZMA_FORMAT_ERROR:
    PANIC("lzma: invalid .xz format (err: %d)", retnum);
  case LZMA_OPTIONS_ERROR:
    PANIC("lzma: invalid options (err: %d)", retnum);
  case LZMA_DATA_ERROR:
  case LZMA_BUF_ERROR:
    PANIC("lzma: file is corrupted or truncated (err: %d)", retnum);
  case LZMA_UNSUPPORTED_CHECK:
    PANIC("lzma: file intergrity check not supported (err: %d)", retnum);
  case LZMA_PROG_ERROR:
    PANIC("lzma: this is bug (err: %d)", retnum);
  }
}

void XzFileHandler::start_lzma_decoder_stream(lzma_stream* stream) {
  constexpr uint64_t memory_limit = std::numeric_limits<uint64_t>::max();
  check_lzma_ret(lzma_stream_decoder(stream, memory_limit, LZMA_TELL_UNSUPPORTED_CHECK | LZMA_CONCATENATED));
}

XzFileHandler::XzFileHandler(const std::string& path, FileMode mode)
    : TextFileHandler(path), m_xz_buffer(LZMA_INTERNAL_BUF_SIZE) {
  const char* openmode;
  switch (mode) {

  case FileMode::READ:
    openmode = "rb";
    start_lzma_decoder_stream(&m_xz_stream);
    break;

  case FileMode::APPEND:
  case FileMode::WRITE:
    PANIC("Only reading is suported.");
    break;
  }

  m_fptr = std::fopen(m_path.c_str(), openmode);
  if (m_fptr == nullptr) {
    lzma_end(&m_xz_stream);
    PANIC("Unable to open file: ", m_path.c_str());
  }
}

XzFileHandler::~XzFileHandler() {
  lzma_end(&m_xz_stream);
  if (m_fptr != nullptr)
    std::fclose(m_fptr);
}

void XzFileHandler::clear() noexcept { std::clearerr(m_fptr); }

void XzFileHandler::seek(uint64_t cursor) {
  // lzma is a stream based compression format, so random access to file is not supported.
  // inefficient implementation: re-decompressing the file from the begining
  // not really an issue if we read from the begining of the file.

  lzma_end(&m_xz_stream);
  m_xz_stream = LZMA_STREAM_INIT;
  start_lzma_decoder_stream(&m_xz_stream);

  // set position to 0
  std::fseek(m_fptr, 0, SEEK_SET);
  constexpr size_t bufsize = LZMA_STREAM_BUF_SIZE;
  char buffer[bufsize];

  while (cursor > bufsize) {
    size_t read_size = this->read(buffer, bufsize);
    assert(read_size == bufsize);
    cursor -= read_size;
  }

  [[maybe_unused]] size_t read_size = this->read(buffer, static_cast<size_t>(cursor));
  assert(read_size == cursor);
}

size_t XzFileHandler::read(char* buffer, size_t size) {

  m_xz_stream.next_out = reinterpret_cast<uint8_t*>(buffer);
  m_xz_stream.avail_out = size;

  while (m_xz_stream.avail_out > 0) {
    if (m_xz_stream.avail_in == 0) {

      m_xz_stream.next_in = m_xz_buffer.data();
      m_xz_stream.avail_in = std::fread(m_xz_buffer.data(), 1, m_xz_buffer.size(), m_fptr);

      if (std::ferror(m_fptr)) {
        PANIC("Something went wrong while reading xz file");
      }
    }

    lzma_action action = std::feof(m_fptr) && m_xz_stream.avail_in == 0 ? LZMA_FINISH : LZMA_RUN;
    lzma_ret ret = lzma_code(&m_xz_stream, action);

    if (ret == LZMA_STREAM_END) {
      break;
    }

    check_lzma_ret(ret);
  }

  return size - m_xz_stream.avail_out;
}

#endif

#ifdef USE_BZIP2

unsigned Bzip2FileHandler::safe_cast(uint64_t size) {
  constexpr size_t max = std::numeric_limits<size_t>::max();
  if (size > max) {
    PANIC("%ld is too big for unsigned in call to bzlib function", size);
  }
  return static_cast<unsigned>(size);
}

void Bzip2FileHandler::check_bz2_retcode(int code) {

  switch (code) {
  case BZ_OK:
  case BZ_RUN_OK:
  case BZ_FLUSH_OK:
  case BZ_FINISH_OK:
  case BZ_STREAM_END:
    return;
  case BZ_SEQUENCE_ERROR:
  case BZ_PARAM_ERROR:
    PANIC("bzip2: bad call to bzlib (err: %d)", code);
  case BZ_MEM_ERROR:
    PANIC("bzip2: memory allocation failed (err: %d)", code);
  case BZ_DATA_ERROR:
    PANIC("bzip2: corrupted file (err: %d)", code);
  case BZ_DATA_ERROR_MAGIC:
    PANIC("bzip2: this file do not seems to be a bz2 file (err: %d)", code);
  // These errors should not occur when using the stream API
  case BZ_CONFIG_ERROR:
    PANIC("bzip2: mis-compiled bzlib (err: %d)", code);
  case BZ_IO_ERROR:
  case BZ_UNEXPECTED_EOF:
  case BZ_OUTBUFF_FULL:
    PANIC("bzip2: unexpected error from bzlib (err: %d)", code);
  default:
    PANIC("bzip2: ???? (err: %d)", code);
  }
}

Bzip2FileHandler::Bzip2FileHandler(const std::string& path, FileMode mode)
    : TextFileHandler(path), m_bz2_buffer(BZ2_INTERNAL_BUF_SIZE) {

  const char* openmode = nullptr;
  switch (mode) {
  case FileMode::READ:
    openmode = "rb";
    m_end_bz2_stream = BZ2_bzDecompressEnd;
    std::memset(&m_bz2_stream, 0, sizeof(bz_stream));
    check_bz2_retcode(BZ2_bzDecompressInit(&m_bz2_stream, 0, 0));
    break;

  case FileMode::APPEND:
  case FileMode::WRITE:
    PANIC("Only reading is suported with bz2 compression format");
    break;
  }

  m_fptr = std::fopen(m_path.c_str(), openmode);
  if (m_fptr == nullptr) {
    m_end_bz2_stream(&m_bz2_stream);
    PANIC("Unable to open file: ", m_path.c_str());
  }
}

Bzip2FileHandler::~Bzip2FileHandler() {
  m_end_bz2_stream(&m_bz2_stream);
  if (m_fptr != nullptr)
    std::fclose(m_fptr);
}


void Bzip2FileHandler::clear() noexcept { std::clearerr(m_fptr); }

void Bzip2FileHandler::seek(uint64_t cursor) {
  // bzip2 is a stream based compression format, so random access to file is not supported.
  // inefficient implementation: re-decompressing the file from the begining
  // not really an issue if we read from the begining of the file.

  m_end_bz2_stream(&m_bz2_stream);
  std::memset(&m_bz2_stream, 0, sizeof(bz_stream));
  check_bz2_retcode(BZ2_bzDecompressInit(&m_bz2_stream, 0, 0));

  // set position to 0
  std::fseek(m_fptr, 0, SEEK_SET);
  constexpr size_t bufsize = BZ2_STREAM_BUF_SIZE;
  char buffer[bufsize];

  while (cursor > bufsize) {
    size_t read_size = read(buffer, bufsize);
    assert(read_size == bufsize);
    cursor -= read_size;
  }

  [[maybe_unused]] size_t read_size = read(buffer, static_cast<size_t>(cursor));
  assert(read_size == cursor);
}

size_t Bzip2FileHandler::read(char* buffer, size_t size) {

  m_bz2_stream.next_out = buffer;
  m_bz2_stream.avail_out = safe_cast(size);

  while (m_bz2_stream.avail_out > 0) {
    if (m_bz2_stream.avail_in == 0) {

      m_bz2_stream.next_in = m_bz2_buffer.data();
      m_bz2_stream.avail_in = safe_cast(std::fread(m_bz2_buffer.data(), 1, m_bz2_buffer.size(), m_fptr));

      if (std::ferror(m_fptr)) {
        PANIC("Something went wrong while reading xz file");
      }
    }

    int retcode = BZ2_bzDecompress(&m_bz2_stream);

    if (retcode == BZ_STREAM_END) {
      break;
    }
    // check for error
    check_bz2_retcode(retcode);
  }

  return size - m_bz2_stream.avail_out;
}

#endif

/* ------------------------------------------------------------------------- */

TextFile::TextFile(std::string filepath, FileMode mode, FileCompression compression)
    : File(std::move(filepath), mode, compression), m_buf(FILE_BUFFER_SIZE, 0), m_lstart(m_buf.data()),
      m_bend(m_buf.data() + m_buf.size()) {

  switch (compression) {
  case FileCompression::BZIP2:
#ifdef USE_BZIP2
    m_handler = std::make_unique<Bzip2FileHandler>(m_path, m_mode);
#endif
    break;

  case FileCompression::XZ:
#ifdef USE_LZMA
    m_handler = std::make_unique<XzFileHandler>(m_path, m_mode);
#endif
    break;

  case FileCompression::GZIP:
#ifdef USE_ZLIB
    m_handler = std::make_unique<GzipFileHandler>(m_path, m_mode);
#endif
    break;

  case FileCompression::NONE:
    m_handler = std::make_unique<ASCIIFileHandler>(m_path, m_mode);
    break;
  }

  if (m_handler == nullptr) {
    PANIC("Unable to instantiate a file handler.");
  }
}

uint64_t TextFile::tell() const {
  uint64_t delta = is_buffer_init() ? static_cast<uint64_t>(m_lstart - m_buf.data()) : 0;
  return m_cursor + delta;
}

void TextFile::seek(uint64_t pos) {
  m_eof = false;
  m_hdlr_eof = false;

  if (is_buffer_init()) {
    int64_t delta = static_cast<int64_t>(pos) - static_cast<int64_t>(m_cursor);
    if (delta >= 0 && delta < static_cast<int64_t>(m_buf.size())) {
      m_lstart = m_buf.data() + delta;
      m_eof = false;
      return;
    }
  }

  m_handler->seek(pos);
  m_cursor = pos;
  m_buf[0] = NULL_CHAR;
}

void TextFile::clear() {
  m_eof = false;
  m_hdlr_eof = false;
  m_handler->clear();
}

void TextFile::reset() {
  clear();
  seek(0);
}

std::string_view TextFile::get_line() {

  // if the buffer is not initialized fill it at pos = 0
  if (!is_buffer_init())
    fill_buffer(0);

  if (m_eof)
    return "";

  size_t length = 0;
  size_t msvc = 0; // windows compatibility

  while (true) {

    // number of charater still in the buffer
    size_t remain = static_cast<size_t>(m_bend - m_lstart);

    // needle is the position of the next new line
    // const char* newline = reinterpret_cast<const char*>(std::memchr(m_lstart + length, '\n', remain - length));
    const char* newline = static_cast<const char*>(std::memchr(m_lstart + length, '\n', remain - length));

    if (newline != nullptr) {
      // if (m_lstart < newline) PANIC("Fail to get_line: lsart < buf.data().");
      length += static_cast<size_t>(newline - m_lstart + 1);

      // check for windows string compatibility
      if (newline > m_buf.data() && newline[-1] == '\r')
        msvc = 1;

      break;

    } else if (m_hdlr_eof) {
      // we reach the end of the file
      m_eof = true;

      if (m_lstart != m_bend - 1) {
        std::string_view line = std::string_view(m_lstart);
        m_lstart += line.length();
        return line;
      }
    }

    if (remain >= m_buf.size()) {
      size_t shift = static_cast<size_t>(m_lstart - m_buf.data());
      m_buf.resize(2 * m_buf.size(), 0);
      m_lstart = m_buf.data() + shift;
      m_bend = m_buf.data() + m_buf.size();
    }

    std::memmove(m_buf.data(), m_lstart, remain);
    fill_buffer(remain);
  }

  std::string_view line = std::string_view(m_lstart, length - msvc - 1);
  m_lstart += length;

  return line;
}

void TextFile::fill_buffer(size_t pos) {
  size_t size = m_buf.size() - pos;

  if (is_buffer_init())
    m_cursor += size;

  size_t read_size = m_handler->read(m_buf.data() + pos, size);

  if (read_size < size) {
    m_hdlr_eof = true;
    std::memset(m_buf.data() + pos + read_size, 0, size - read_size);
  }

  m_lstart = m_buf.data();
}

/* ------------------------------------------------------------------------- */

void TextParser::scan() {
  // if end of file, return.
  if (m_eof)
    return;

  size_t i = m_file.tell();

  while (!m_file.eof()) {
    if (int64_t c = next(); c > -1) {
      m_loc.push_back(static_cast<size_t>(c));
    }
    break;
  }

  m_eof = true;
  m_file.clear();

  if (i == 0 && !m_loc.empty()) {
    m_file.seek(m_loc[0]);
  } else {
    m_file.seek(i);
  }
}

bool TextParser::operator()(IOContext& ctx) {
  m_file.seek(0);
  return parse(ctx);
}

bool TextParser::operator()(IOContext& ctx, size_t index) {

  if (index >= size())
    scan();

  if (size() == 0) {
    return false;
  }

  if (index >= size()) {
    return false;
  }

  m_file.seek(static_cast<uint64_t>(m_loc[index]));
  return parse(ctx);
}

/* ------------------------------------------------------------------------- */

ParserFactory::ParserFactory() {
  register_format<lmpdata::LAMMPSDataParser>();
  register_format<lmpdump::LAMMPSDumpParser>();
  register_format<xyz::ExtendedXYZParser>();
};

std::string ParserFactory::infos() const {

  std::string infos = "Supported atom format:\n";

  std::unordered_map<size_t, RegisteredFormat>::const_iterator it;
  size_t size = index_to_parser.size();

  for (size_t index = 0; index < size; ++index) {
    infos += strf("\t- %s\n", index_to_parser.find(index)->second.metadata.name.c_str());
  }

  infos += "\n";
  return infos;
};

const ParserFactory& parser_factory() {
  static const ParserFactory parser_factory_instance;
  return parser_factory_instance;
};

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lmpdata {

int64_t LAMMPSDataParser::next() {
  // LAMMPS Data only contains one frame
  if (size_t cursor = m_file.tell(); cursor == 0) {
    m_file.get_line();
    return static_cast<int64_t>(cursor);
  } else {
    return -1;
  }
}

bool LAMMPSDataParser::parse(IOContext& ctx) {

  while (!m_file.eof()) {
    if (!skipline(current_line())) {

      switch (m_current_section) {
      case Section::Header:
        parse_section_header(ctx);
        if (ctx.has<IOContext::DOMAIN_ONLY>()) {
          return true;
        }
        break;
      case Section::Masses:
        parse_section_masses();
        break;
      case Section::Atoms:
        parse_section_atoms(ctx);
        break;
      case Section::Velocities:
        parse_section_velocities(ctx);
        break;
      case Section::Ignored:
        jumpto_next_section();
        break;
      default:
        break;
      }
    }

    if (m_current_section != Section::Ignored) {
      current_line() = trim_spaces(m_file.get_line());
    }
  }
  return true;
}

bool LAMMPSDataParser::parse_section_header(IOContext& ctx) {

  std::array<bool, enum_size<Keyword>()> matched_keywords{};

  while (!m_file.eof()) {

    if (!skipline(current_line())) {

      tokenize(current_line(), m_tokens);

      for (size_t i = 0; i < enum_size<Keyword>(); ++i) {
        if (matched_keywords[i]) {
          continue;
        }

        if ((*keyword_parsers[i])(m_tokens, ctx)) {
          matched_keywords[i] = true;
          break;
        }
      }

      // if any section is parsed we break from parsing header
      if (matched_keywords[enum_to_index(Keyword::Section)]) {
        break;
      }
    }

    current_line() = trim_spaces(m_file.get_line());
  }

  m_current_section = Section::Ignored;
  return true;
}

bool LAMMPSDataParser::parse_section_atoms(IOContext& ctx) {

  size_t index, style;
  tokenize(section_line(), m_tokens);
  if (token_set_match_any(m_tokens, NEEDLES_ATOM_STYLE, index, style)) {
    m_atom_style = enum_from_index<AtomStyle>(style);
  } else {
    PANIC("Invalid atom style defined.");
  }

  size_t particle_count = 0;
  ctx.data.positions.resize(ctx.data.particle_count);
  ctx.data.velocities.resize(ctx.data.particle_count);
  auto& parser =
      (ctx.has<IOContext::REMAP_ATOM>() ? atom_data_parsers<true>() : atom_data_parsers<false>())[style];

  while (!m_file.eof()) {

    if (!skipline(current_line())) {
      tokenize(current_line(), m_tokens);

      if (parse_keyword_section(m_tokens, ctx)) {
        break;
      }

      bool atom_overflow = cmp::gt(++particle_count, ctx.data.particle_count);
      if (atom_overflow || !parser(m_tokens, ctx, particle_count)) {
        if (atom_overflow) {
          PANIC("Too many atoms expected %ld got %ld", ctx.data.particle_count, particle_count);
        }
        PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
      }
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  if (particle_count != ctx.data.particle_count) {
    PANIC("Number of particles is not consitent...");
  }

  m_current_section = Section::Ignored;
  return true;
}

// Currently we ignore the mass section
bool LAMMPSDataParser::parse_section_masses() {
  m_current_section = Section::Ignored;
  return true;
}

bool LAMMPSDataParser::parse_section_velocities(IOContext& ctx) {

  size_t particle_count = 0;
  auto parser = ctx.has<IOContext::REMAP_ATOM>() ? parse_atom_velocities<true> : parse_atom_velocities<false>;

  while (!m_file.eof()) {

    if (!skipline(current_line())) {
      tokenize(current_line(), m_tokens);

      if (parse_keyword_section(m_tokens, ctx)) {
        break;
      }

      bool atom_overflow = cmp::gt(++particle_count, ctx.data.particle_count);
      if (atom_overflow || !parser(m_tokens, ctx, particle_count)) {
        if (atom_overflow) {
          PANIC("Too many atoms expected %ld got %ld", ctx.data.particle_count, particle_count);
        }
        PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
      }
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  if (particle_count != ctx.data.particle_count) {
    PANIC("Number of particles is not consitent...");
  }

  m_current_section = Section::Ignored;
  return true;
}

void LAMMPSDataParser::jumpto_next_section() {
  current_line() = trim_spaces(current_line());
  while (!m_file.eof()) {
    if (!skipline(current_line())) {

      tokenize(current_line(), m_tokens);

      if (token_match_any(m_tokens[0], NEEDLES_SECTIONS)) {
        set_current_section(m_tokens[0]);
        section_line() = current_line(); // keep track of the section for later
        break;
      }
    }
    current_line() = trim_spaces(m_file.get_line());
  }
}

bool parse_keyword_atoms(const TokenSet& tokens, IOContext& ctx) {
  return tokens.atleast(2) && token_match_any(tokens.at<1>(), NEEDLES_KEYWORD_ATOMS) &&
         parse_tokens(tokens.at<0>(), ctx.data.particle_count);
}

bool parse_keyword_atom_types(const TokenSet& tokens, IOContext& ctx) {
  return tokens.atleast(3) && token_set_match_seq(tokens, NEEDLES_KEYWORD_ATOM_TYPES) &&
         parse_tokens(tokens.at<0>(), ctx.data.type_count);
}

template <size_t I, size_t J>
bool parse_keyword_box_axis(const TokenSet& tokens, IOContext& ctx, const TokenNeedles<2>& needles) {
  double boxlo, boxhi;
  if (tokens.atleast(4) && token_set_match_seq(tokens, needles) &&
      parse_tokens(tokens.at<0>(), boxlo, tokens.at<1>(), boxhi)) {
    matrix_traits::at<I, J>(ctx.data.cell) = boxhi - boxlo;
    return true;
  }
  return false;
}

inline bool parse_keyword_xaxis(const TokenSet& tokens, IOContext& ctx) {
  return parse_keyword_box_axis<0, 0>(tokens, ctx, NEEDLES_KEYWORD_XAXIS);
}

inline bool parse_keyword_yaxis(const TokenSet& tokens, IOContext& ctx) {
  return parse_keyword_box_axis<1, 1>(tokens, ctx, NEEDLES_KEYWORD_YAXIS);
}

inline bool parse_keyword_zaxis(const TokenSet& tokens, IOContext& ctx) {
  return parse_keyword_box_axis<2, 2>(tokens, ctx, NEEDLES_KEYWORD_ZAXIS);
}

bool parse_keyword_tilts(const TokenSet& tokens, IOContext& ctx) {
  double* xy = &ctx.data.cell.m12;
  double* xz = &ctx.data.cell.m13;
  double* yz = &ctx.data.cell.m23;
  return tokens.atleast(6) && token_set_match_seq(tokens, NEEDLES_KEYWORD_TILTS) &&
         parse_tokens(tokens.at<0>(), *xy, tokens.at<1>(), *xz, tokens.at<2>(), *yz);
}

inline bool parse_keyword_section(const TokenSet& tokens, IOContext&) {
  return tokens.atleast(1) && token_match_any(tokens.at<0>(), NEEDLES_SECTIONS);
}

} // namespace lmpdata


/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Dump format

namespace lmpdump {

int64_t LAMMPSDumpParser::next() {
  if (size_t cursor = m_file.tell(); cursor == 0) {
    m_file.get_line();
    return static_cast<int64_t>(cursor);
  } else {
    return -1;
  }
}

bool LAMMPSDumpParser::parse(IOContext& ctx) {

  constexpr StringView label_natom = "ITEM: NUMBER OF ATOMS";
  constexpr StringView label_box_bounds = "ITEM: BOX BOUNDS";
  constexpr StringView label_atoms = "ITEM: ATOMS";

  Properties properties;


  while (!m_file.eof() && (skipline(current_line()) || !startswith(current_line(), label_natom))) {
    current_line() = trim_spaces(m_file.get_line());
  }


  current_line() = trim_spaces(m_file.get_line());
  tokenize(current_line(), m_tokens, 1);
  if (!parse_tokens(m_tokens.at<0>(), ctx.data.particle_count)) {
    PANIC("Invalid number of atom: %s", std::string(current_line()).c_str());
  }


  // 'ITEM: BOX BOUNDS' is expected just after NUMBER OF ATOMS
  current_line() = trim_spaces(m_file.get_line());
  if (!startswith(current_line(), label_box_bounds)) {
    PANIC("");
  }


  tokenize(current_line(), m_tokens);
  properties.triclinic = m_tokens.atleast(7) ? true : false;
  Mat3d bbox;
  for (size_t i = 0; i < 3; ++i) {
    tokenize(trim_spaces(m_file.get_line()), m_tokens);
    if (!(properties.triclinic
              ? parse_box_bounds_axis<true>(m_tokens, matrix_traits::row(i, bbox))
              : parse_box_bounds_axis<false>(m_tokens, matrix_traits::row(i, bbox)))) {
      PANIC("");
    }
  }
  convert_bbox_to_lattice_vectors(ctx, bbox);

  if (ctx.has<IOContext::DOMAIN_ONLY>()) {
    return true;
  }

  // 'ITEM: ATOMS' is expected just after BOX BOUNDS
  current_line() = trim_spaces(m_file.get_line());
  if (!startswith(current_line(), label_atoms)) {
    PANIC("");
  }

  tokenize(current_line(), m_tokens);
  (properties.triclinic)
      ? parse_atom_properties<true>(m_tokens, properties)
      : parse_atom_properties<false>(m_tokens, properties);

  // parse atom data
  size_t particle_count = 0;
  ctx.data.positions.resize(ctx.data.particle_count);
  ctx.data.velocities.resize(ctx.data.particle_count);
  auto parser =
      (ctx.has<IOContext::REMAP_ATOM>() && properties.has<Field::ID>())
          ? parse_atom_data<true>
          : parse_atom_data<false>;

  current_line() = trim_spaces(m_file.get_line());

  while (!m_file.eof() && particle_count < ctx.data.particle_count) {

    if (!skipline(current_line())) {
      tokenize(current_line(), m_tokens);

      bool atom_overflow = cmp::gt(++particle_count, ctx.data.particle_count);
      if (atom_overflow || !parser(m_tokens, ctx, properties, particle_count)) {
        if (atom_overflow) {
          PANIC("Too many atoms expected %ld got %ld", ctx.data.particle_count, particle_count);
        }
        PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
      }
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  if (particle_count != ctx.data.particle_count) {
    PANIC("Number of particles is not consitent...");
  }

  return true;
}

} // namespace lmpdump

/* ------------------------------------------------------------------------- */

namespace xyz {

int64_t ExtendedXYZParser::next() {
  if (size_t cursor = m_file.tell(); cursor == 0) {
    m_file.get_line();
    return static_cast<int64_t>(cursor);
  } else {
    return -1;
  }
}

bool ExtendedXYZParser::parse(IOContext& ctx) {

  while (!m_file.eof() && skipline(current_line())) {
    current_line() = trim_spaces(m_file.get_line());
  }

  // check if the first non empty line is a valid number
  tokenize(current_line(), m_tokens);
  if (!parse_tokens(m_tokens.at<0>(), ctx.data.particle_count)) {
    PANIC("Invalid atom number: '%s'", std::string(current_line()).c_str());
  }

  // next line should be the comment line with lattice and properties definitions
  current_line() = trim_spaces(m_file.get_line());

  if (!parse_xyz_comment_line(current_line(), m_tokens, m_ptokens) ||
      !(parse_xyz_lattice(m_tokens, ctx) && parse_xyz_properties(m_ptokens, m_properties))) {
    PANIC("Failed to read the comment line: '%s'", std::string(current_line()).c_str());
  }

  if (ctx.has<IOContext::DOMAIN_ONLY>()) {
    return true;
  }

  // parse atom data
  ctx.data.positions.resize(ctx.data.particle_count);
  ctx.data.velocities.resize(ctx.data.particle_count);
  size_t particle_count = 0;

  auto parser =
      (ctx.has<IOContext::REMAP_ATOM>() && m_properties.has<Field::ID>())
          ? parse_xyz_atom_line<true>
          : parse_xyz_atom_line<false>;

  current_line() = trim_spaces(m_file.get_line());

  while (!m_file.eof() && particle_count < ctx.data.particle_count) {

    if (!skipline(current_line())) {
      tokenize(current_line(), m_tokens);

      bool atom_overflow = cmp::gt(++particle_count, ctx.data.particle_count);
      if (atom_overflow || !parser(m_tokens, ctx, m_properties, particle_count)) {
        if (atom_overflow) {
          PANIC("Too many atoms expected %ld got %ld", ctx.data.particle_count, particle_count);
        }
        PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
      }
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  if (particle_count != ctx.data.particle_count) {
    PANIC("Number of particles is not consitent...");
  }

  return true;
}
} // namespace xyz

/* ------------------------------------------------------------------------- */

std::string remove_leading_chars(const std::string& str, const char c) {
  size_t start = 0;
  while (start < str.size() && str[start] == c) {
    ++start;
  }
  return str.substr(start);
}

const FileInfos guess_file_infos(const std::string& filepath) {

  std::string extension;
  std::filesystem::path path(filepath);

  // remove leading dots return by extension()
  extension = remove_leading_chars(path.extension(), '.');

  // if no extension can't do anything here
  if (extension.empty()) {
    return {.compression = FileCompression::NONE, .format = parser_factory().from_alias("")};
  }

  FileCompression compression = convert_char_to_compression(extension);

  // if first extension is related to compression, get the second one
  if (compression != FileCompression::NONE) {
    extension = remove_leading_chars(path.stem().extension(), '.');
  }

  // return file infos either if format is invalid, this is handled after
  return {.compression = compression, .format = parser_factory().from_alias(extension)};
}

const FileInfos get_file_infos(const std::string& filepath, std::string user_format, std::string user_compression) {

  // never trust the user input
  FileInfos guessed_infos = guess_file_infos(filepath);
  RegisteredFormat format = parser_factory().from_alias(user_format);

  if ((!format.is_valid) && (!guessed_infos.format.is_valid)) {
    PANIC("Could not guess file format.");
  }

  FileInfos infos{.compression = user_compression.empty()
                                     ? guessed_infos.compression
                                     : convert_char_to_compression(user_compression),
                  .format = format.is_valid ? format : guessed_infos.format};
  return infos;
}

}
