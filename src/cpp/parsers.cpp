#include <filesystem>
#include "parsers.h"

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
      size_t shift = m_lstart - m_buf.data();
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
      loc.push_back(c);
    }
    break;
  }

  m_eof = true;
  m_file.clear();

  if (i == 0 && !loc.empty()) {
    m_file.seek(loc[0]);
  } else {
    m_file.seek(i);
  }
}

bool TextParser::operator()(Context& ctx) {
  m_file.seek(0);
  return parse(ctx);
}

bool TextParser::operator()(Context& ctx, size_t index) {

  if (index >= size())
    scan();

  if (size() == 0) {
    return false;
  }

  if (index >= size()) {
    return false;
  }

  m_file.seek(static_cast<uint64_t>(loc[index]));
  return parse(ctx);
}

/* ------------------------------------------------------------------------- */

ParserFactory::ParserFactory() {
  register_format<lmpdata::LAMMPSDataParser>();
  // register_format<LAMMPSDumpParser>();
  // register_format<xyz::ExtendedXYZParser>();
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

// LAMMPS Data only contains one frame
int64_t LAMMPSDataParser::next() {
  if (size_t cursor = m_file.tell(); cursor == 0) {
    m_file.get_line();
    return static_cast<int64_t>(cursor);
  } else {
    return -1;
  }
}


bool LAMMPSDataParser::parse(Context& ctx) {

  // lout << "Parsing LAMMPS Data file (v2):" << std::endl;

  current_line() = m_file.get_line();

  while (!m_file.eof()) {

    current_line() = trim_spaces(current_line());
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }

    switch (m_current_section) {
    case Section::HEADER:
      parse_header(ctx);
      // we don't care about the rest of the file
      return true;
      break;
    case Section::MASSES:
      parse_masses();
      break;
    case Section::ATOMS:
      parse_atoms(ctx);
      break;
    case Section::IGNORED:
      forward_to_next_section();
      break;
    default:
      break;
    }

    current_line() = m_file.get_line();
  }

  return true;
}

bool LAMMPSDataParser::parse_header(Context& ctx) {

  std::array<bool, static_cast<size_t>(Field::END)> matched_field = {};

  while (!m_file.eof()) {

    current_line() = trim_spaces(current_line());
    if (skipline(current_line())) {
      current_line() = trim_spaces(m_file.get_line());
      continue;
    }
    tokenize(current_line(), m_tokens);

    for (size_t i = 0; i < field_readers().size(); ++i) {

      if (matched_field[i])
        continue;

      // auto& [field, parser] = m_field_parsers[i];
      auto& [field, parser] = field_readers()[i];

      if (parser(m_tokens, ctx)) {
        matched_field[i] = true;
        break;
      }
    }

    // if any section keyword is meet we break from parsing headers
    // and go to the next section
    if (matched_field[6]) {
      m_current_section = get_section(m_tokens[0]);
      break;
    }
    current_line() = m_file.get_line();
  }
  return true;
}

bool LAMMPSDataParser::parse_atoms(Context& ctx) {
  m_current_section = Section::IGNORED;
  return true;
}

bool LAMMPSDataParser::parse_masses() {
  m_current_section = Section::IGNORED;
  return true;
}

bool LAMMPSDataParser::parse_velocities(Context&) {
  m_current_section = Section::IGNORED;
  return true;
}

void LAMMPSDataParser::forward_to_next_section() {
  while (!m_file.eof()) {

    current_line() = trim_spaces(current_line());
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }

    tokenize(current_line(), m_tokens);

    if (token_match_any(m_tokens[0], tok_sections)) {
      m_current_section = get_section(m_tokens[0]);
      section_line() = current_line();
      break;
    }

    current_line() = m_file.get_line();
  }
}

bool read_field_natom(const TokenSet& tokens, Context& ctx) {
  return token_match_any(tokens[1], tok_atoms) && (tokens.size() >= 2) &&
         tokens_to_num(tokens[0], ctx.n_particles);
}

bool read_field_atom_types(const TokenSet& tokens, Context& ctx) {
  if (!token_set_match_seq(tokens, tok_atom_types))
    return false;
  return tokens_to_num(tokens[0], ctx.n_species);
}

bool read_field_xaxis(const TokenSet& tokens, Context& ctx) {
  double xlo, xhi;
  bool ok = token_set_match_seq(tokens, tok_xaxis) && tokens_to_num(tokens[0], xlo, tokens[1], xhi);
  if (ok) {
    ctx.cell[0][0] = xhi - xlo;
  }
  return ok;
}

bool read_field_yaxis(const TokenSet& tokens, Context& ctx) {
  double ylo, yhi;
  bool ok = token_set_match_seq(tokens, tok_yaxis) && tokens_to_num(tokens[0], ylo, tokens[1], yhi);
  if (ok) {
    ctx.cell[1][1] = yhi - ylo;
  }
  return ok;
}

bool read_field_zaxis(const TokenSet& tokens, Context& ctx) {
  double zlo, zhi;
  bool ok = token_set_match_seq(tokens, tok_zaxis) && tokens_to_num(tokens[0], zlo, tokens[1], zhi);
  if (ok) {
    ctx.cell[2][2] = zhi - zlo;
  }
  return ok;
}

bool read_field_tilt(const TokenSet& tokens, Context& ctx) {
  double xy, xz, yz;
  bool ok = token_set_match_seq(tokens, tok_tilt) && (tokens.size() >= 6) &&
            tokens_to_num(tokens[0], xy, tokens[1], xz, tokens[2], yz);
  if (ok) {
    ctx.cell[0][1] = xy;
    ctx.cell[0][2] = xz;
    ctx.cell[1][2] = yz;
  }
  return ok;
}

bool read_field_section(const TokenSet& tokens, Context& ctx) { 
  return token_set_match_any(tokens, tok_sections);
}

} // namespace lmpdata

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
    PANIC("AAAAAAAAA");
  }

  FileInfos infos{.compression = user_compression.empty()
                                     ? guessed_infos.compression
                                     : convert_char_to_compression(user_compression),
                  .format = format.is_valid ? format : guessed_infos.format};
  return infos;
}

}
