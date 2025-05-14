/* ----------------------------------------------------------------------------
Parser function to read atomic configuration file. Light weight parser copied
from exaStamp parser to only parse lattice vector.
---------------------------------------------------------------------------- */
#pragma once

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <functional>

#include "tokens.h"

#define FILE_BUFFER_SIZE 8192
#define NULL_CHAR '\0'
#define LZMA_STREAM_BUF_SIZE 8192
#define LZMA_INTERNAL_BUF_SIZE 8192
#define BZ2_STREAM_BUF_SIZE 8192
#define BZ2_INTERNAL_BUF_SIZE 8192

#define PANIC(msg) throw std::runtime_error(msg);

namespace io {
using namespace tokens;

/* ------------------------------------------------------------------------- */

// Struct that hold lattice vectors
struct Context {
  size_t n_particles = 0;
  size_t n_species = 0;
  double cell[3][3];
};

/* ------------------------------------------------------------------------- */

// File openning mode
enum FileMode : char {
  READ = 'r',
  WRITE = 'w',
  APPEND = 'a',
};

// File compression mode
enum FileCompression {
  NONE,
  GZIP,
  BZIP2,
  XZ,
};

/* ------------------------------------------------------------------------- */

inline const char* convert_mode_to_char(FileMode mode) {
  switch (mode) {
  case FileMode::READ:
    return "rb";
    break;
  case FileMode::APPEND:
    return "a+b";
    break;
  case FileMode::WRITE:
    return "wb";
    break;
  default:
    PANIC("Invalid convertion from Filemode to char *");
    break;
  }
}

inline FileCompression convert_char_to_compression(const std::string& str) {
  if (str == "gz")
    return FileCompression::GZIP;
  if (str == "bz2")
    return FileCompression::BZIP2;
  if (str == "xz")
    return FileCompression::XZ;
  return FileCompression::NONE;
}

inline std::string convert_compression_to_char(FileCompression& compression) {
  if (compression == FileCompression::GZIP)
    return "gz";
  if (compression == FileCompression::BZIP2)
    return "bz2";
  if (compression == FileCompression::XZ)
    return "xz";
  return "";
}

/* ------------------------------------------------------------------------- */

// Base class that handler basic file operation
// Specialization allow abstraction on file types.
class TextFileHandler {
public:
  TextFileHandler(const std::string& filepath);

  virtual ~TextFileHandler() = default;

  virtual void clear() noexcept = 0;
  virtual void seek(uint64_t position) = 0;
  virtual size_t read(char* buffer, size_t size) = 0;

protected:
  const std::string& path() { return m_path; };
  std::string m_path;
};

class ASCIIFileHandler final : public TextFileHandler {
public:
  ASCIIFileHandler(const std::string& path, FileMode mode);
  ~ASCIIFileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* data, size_t size) override;

private:
  std::FILE* m_fptr;
};

#ifdef USE_ZLIB
#define ZLIB_CONST
#include <zconf.h>
#include <zlib.h>

typedef struct gzFile_s* gzFile;

class GzipFileHandler final : public TextFileHandler {
public:
  GzipFileHandler(const std::string& path, FileMode mode);
  ~GzipFileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* buffer, size_t size) override;

  const char* gz_error() const;
  static unsigned safe_cast(size_t);

private:
  gzFile m_fptr = nullptr;
};

#endif

#ifdef USE_LZMA
#include <lzma.h>

class XzFileHandler final : public TextFileHandler {
public:
  XzFileHandler(const std::string& path, FileMode mode);
  ~XzFileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* buffer, size_t size) override;

  size_t safe_cast(uint64_t value);
  void check_lzma_ret(lzma_ret ret);
  void start_lzma_decoder_stream(lzma_stream* stream);

private:
  std::FILE* m_fptr = nullptr;
  lzma_stream m_xz_stream = LZMA_STREAM_INIT;

  std::vector<uint8_t> m_xz_buffer = {0};
  // static constexpr size_t m_xz_buffer_end = IOEXTRA_LZMA_INTERNAL_BUF_SIZE;
  // uint8_t* m_xz_buffer[m_xz_buffer_end];
};

#endif

#ifdef USE_BZIP2
#include <bzlib.h>

class Bzip2FileHandler final : public TextFileHandler {
public:
  Bzip2FileHandler(const std::string& path, FileMode mode);
  ~Bzip2FileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* buffer, size_t size) override;

  unsigned safe_cast(uint64_t size);
  void check_bz2_retcode(int code);

private:
  std::FILE* m_fptr;
  std::function<int(bz_stream*)> m_end_bz2_stream;
  bz_stream m_bz2_stream;
  std::vector<char> m_bz2_buffer;
};

#endif

/* ------------------------------------------------------------------------- */

class File {
protected:
  std::string m_path;
  FileMode m_mode;
  FileCompression m_compression;

  File(std::string path, FileMode mode, FileCompression compression)
      : m_path(std::move(path)), m_mode(mode), m_compression(compression) {}

public:
  inline const std::string& path() const { return m_path; };
  inline FileMode mode() const { return m_mode; };
  inline FileCompression compression() const { return m_compression; };
};

class TextFile : public File {
private:
  std::vector<char> m_buf;
  const char* m_lstart; // line start
  const char* m_bend;   // buf end

  uint64_t m_cursor = 0;
  bool m_eof = false;
  bool m_hdlr_eof = false;

  std::unique_ptr<TextFileHandler> m_handler;

  inline bool is_buffer_init() const { return m_buf[0] != NULL_CHAR; }
  void fill_buffer(size_t pos);

public:
  TextFile(std::string filepath, FileMode mode, FileCompression compression);

  inline bool eof() const { return m_eof; };

  uint64_t tell() const;
  void seek(uint64_t pos);
  void clear();
  void reset();

  // get the next line in the buffer.
  // If the line is not complete, load the next section of the buffer
  std::string_view get_line();
};

/* ------------------------------------------------------------------------- */

struct Metadata {
  const std::string name = "";
  const std::vector<std::string> aliases = {};
  const std::string description = "";
};

template<typename T>
const Metadata get_parser_metadata() { 
  return { };
}

// Base class to all format parsers
class Parser {
public:
  Parser() = default;
  virtual ~Parser() = default;

  Parser(const Parser&) = delete;
  Parser& operator=(const Parser&) = delete;
  Parser(Parser&&) = delete;
  Parser& operator=(Parser&&) = delete;

  // read the first step in the file
  virtual inline bool operator()(Context&) { return false; };
  // read the step at the given index
  virtual inline bool operator()(Context&, size_t) { return false; };
  // returns the number of known step in the file
  virtual size_t size() = 0;
};

// Specialization for format that are just text files
class TextParser : public Parser {

protected:
  bool m_eof = false;
  std::vector<size_t> loc{};
  TextFile m_file;

  std::array<std::string_view, 10> line_buffer{};
  inline std::string_view& current_line(void) { return line_buffer[0]; };
  TokenSet m_tokens{};

  // std::string_view m_current_line;

public:
  TextParser(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression) {}

  virtual ~TextParser() override = default;

  virtual int64_t next() = 0;
  virtual bool parse(Context& ctx) = 0;

  bool operator()(Context& ctx) override;
  bool operator()(Context& ctx, size_t index) override;

  void scan();

  inline size_t size() override { return loc.size(); };
  inline const std::string_view& current_line(void) const { return line_buffer[0]; }
};

/* ------------------------------------------------------------------------- */

using ParserCreator = std::function<std::unique_ptr<Parser>(std::string, FileMode, FileCompression)>;
// using ParserCreator = std::unique_ptr<Parser>(*)(std::string, FileMode, FileCompression);

struct RegisteredFormat {
  const Metadata metadata;
  const ParserCreator creator;
  bool is_valid = true;
};

class ParserFactory {
public:

  ParserFactory();

  std::unordered_map<std::string, size_t> alias_to_index;
  std::unordered_map<size_t, RegisteredFormat> index_to_parser;

  template <typename T> void register_format() {

    size_t index = index_to_parser.size();
    Metadata metadata = get_parser_metadata<T>();

    index_to_parser.insert(
        {index,
         {.metadata = metadata,
          .creator = [](std::string filepath, FileMode mode, FileCompression compression) {
            return std::make_unique<T>(filepath, mode, compression);
          }}});

    // map extension to the index
    for (const std::string& alias : metadata.aliases) {
      std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
      if (it != alias_to_index.end()) {
        continue; // the alias is already used
      }
      alias_to_index.insert({alias, index});
    }
  }

  RegisteredFormat from_alias(const std::string& alias) const {
    std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
    if (it == alias_to_index.end()) {
      return {.metadata = Metadata{}, .creator = nullptr, .is_valid = false};
    }
    return index_to_parser.find(it->second)->second;
  }

  inline bool has_alias(const std::string& alias) const {
    return alias_to_index.find(alias) != alias_to_index.end();
  }

  std::string infos() const;
};

const ParserFactory& parser_factory();

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format
//
namespace lmpdata {

enum class AtomStyle { ATOMIC, MOLECULAR, FULL, NONE };

enum class Field { NATOM, TYPAT, XAXIS, YAXIS, ZAXIS, TILT, SECTION, END };

enum class Section {
  HEADER,
  MASSES,
  ATOMS,
  VELOCITIES,
  IGNORED,
};

constexpr std::array<Section, 3> sections_list { Section::MASSES, Section::ATOMS, Section::VELOCITIES };

constexpr TokenNeedles<1> tok_atoms = {"atoms"};
constexpr TokenNeedles<2> tok_atom_types = {"atom", "types"};
constexpr TokenNeedles<2> tok_xaxis = {"xlo", "xhi"};
constexpr TokenNeedles<2> tok_yaxis = {"ylo", "yhi"};
constexpr TokenNeedles<2> tok_zaxis = {"zlo", "zhi"};
constexpr TokenNeedles<3> tok_tilt = {"xy", "xz", "yz"};
constexpr TokenNeedles<3> tok_sections = {"Masses", "Atoms", "Velocities"};
constexpr TokenNeedles<3> tok_atom_style = {"atomic", "full", "molecule"};


class LAMMPSDataParser : public TextParser {
public:
  LAMMPSDataParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression), m_atom_style(AtomStyle::NONE),
        m_current_section(Section::HEADER) {}

  int64_t next() override;
  bool parse(Context& proxy) override;

private:

  AtomStyle m_atom_style;
  Section m_current_section;
  
  // keep track of the current section line
  inline std::string_view& section_line(void) { return line_buffer[1]; };
  inline const std::string_view& section_line(void) const { return line_buffer[1]; };

  bool parse_header(Context& ctx);
  bool parse_atoms(Context& ctx);
  bool parse_masses();
  bool parse_velocities(Context& ctx);

  void forward_to_next_section();
};

using AtomLineReader = bool(*)(const TokenSet&, Context&, size_t);
using FieldParser = bool(*)(const TokenSet&, Context&);

bool read_atom_line_atomic(const TokenSet&, Context&, size_t);

inline const std::unordered_map<AtomStyle, AtomLineReader>& atom_line_readers() {
  static std::unordered_map<AtomStyle, AtomLineReader> map = {
    { AtomStyle::ATOMIC , &read_atom_line_atomic }
  };
  return map;
}

bool read_field_natom(const TokenSet&, Context&);
bool read_field_atom_types(const TokenSet&, Context&);
bool read_field_xaxis(const TokenSet&, Context&);
bool read_field_yaxis(const TokenSet&, Context&);
bool read_field_zaxis(const TokenSet&, Context&);
bool read_field_tilt(const TokenSet&, Context&);
bool read_field_section(const TokenSet&, Context&);

inline const std::array<std::pair<Field, FieldParser>, 7>& field_readers() {
  static constexpr std::array<std::pair<Field, FieldParser>, 7> readers = {{
      {Field::NATOM, &read_field_natom},
      {Field::TYPAT, &read_field_atom_types},
      {Field::XAXIS, &read_field_xaxis},
      {Field::YAXIS, &read_field_zaxis},
      {Field::ZAXIS, &read_field_yaxis},
      {Field::TILT, &read_field_tilt},
      {Field::SECTION, &read_field_section},
  }};
  return readers;
}

inline Section get_section(std::string_view str) {
  size_t index;
  if (token_match_any(str, tok_sections, index))
    return sections_list[index];
  else
    return Section::IGNORED;
}

} // namespace lmpdata

template<> inline const Metadata get_parser_metadata<lmpdata::LAMMPSDataParser>() {
  return {
    .name        = "LAMMPS Data",
    .aliases     = { "lmp", "lmp-data" },
    .description = "LAMMPS Data format"
  };
}

/* ------------------------------------------------------------------------- */

struct FileInfos {
  FileCompression compression = FileCompression::NONE;
  const RegisteredFormat format{};
};

std::string remove_leading_chars(const std::string& str, const char c);
const FileInfos guess_file_infos(const std::string& filepath);
const FileInfos get_file_infos(const std::string& filepath, std::string format, std::string compression);

inline Context read_atom_file(const std::string& file, const std::string& format = "", const std::string& compression = "") {

  // Ensure the file exists
  std::string filepath = std::string(std::filesystem::absolute(file));
  if (!std::filesystem::exists(filepath)) {
    PANIC(strf("Input file doest not exists: '%s'\n", filepath.c_str()));
  }

  // create reader object
  FileInfos infos = get_file_infos(filepath, format, compression);
  std::unique_ptr<Parser> parser =
      infos.format.creator(filepath, FileMode::READ, infos.compression);
  if (parser == nullptr) {
    PANIC("Something went wrong... I AM PANICKING");
  }

  Context ctx{};
  if (!(*parser)(ctx)) {
    PANIC("Something went wrong when parsing the file");
  }
  return ctx;
}

} // namespace io
