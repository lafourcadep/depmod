/* ----------------------------------------------------------------------------
Parser function to read atomic configuration file.
---------------------------------------------------------------------------- */
#pragma once

#include "tokens.h"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#define READ_BUFFER_SIZE 1'048'576
#define FILE_BUFFER_SIZE READ_BUFFER_SIZE
#define LZMA_STREAM_BUF_SIZE READ_BUFFER_SIZE
#define LZMA_INTERNAL_BUF_SIZE READ_BUFFER_SIZE
#define BZ2_STREAM_BUF_SIZE READ_BUFFER_SIZE
#define BZ2_INTERNAL_BUF_SIZE READ_BUFFER_SIZE
#define NULL_CHAR '\0'

#define PANIC(...) throw std::runtime_error(strf(__VA_ARGS__));

namespace io {
using namespace tokens;
using namespace tokens::numeric;
using namespace enum_traits;

struct IJK {
  size_t i;
  size_t j;
  size_t k;
};

struct Vec3d {
  double x;
  double y;
  double z;
};

struct Mat3d {
  double m11, m12, m13;
  double m21, m22, m23;
  double m31, m32, m33;
};

namespace matrix_traits {

struct Mat3dSlice {
  double &x, &y, &z;
};

template <size_t N>
inline constexpr Mat3dSlice column(Mat3d& m) {
  static_assert(N >= 0 && N < 2);
  double* base = &m.m11;
  return Mat3dSlice{base[N], base[N + 3], base[N + 6]};
};

template <size_t N>
inline constexpr Mat3dSlice row(Mat3d& m) {
  static_assert(N >= 0 && N < 2);
  double* base = &m.m11;
  size_t index = 3 * N;
  return Mat3dSlice{base[index], base[index + 1], base[index + 2]};
};

inline Mat3dSlice column(size_t i, Mat3d& m) {
  double* base = &m.m11;
  return Mat3dSlice{base[i], base[i + 3], base[i + 6]};
};

inline Mat3dSlice row(size_t i, Mat3d& m) {
  double* base = &m.m11;
  size_t index = 3 * i;
  return Mat3dSlice{base[index], base[index + 1], base[index + 2]};
};

template<size_t I, size_t J>
constexpr inline double& at(Mat3d& m) {
  double* base = &m.m11;
  static_assert(I >= 0 && I < 3);
  static_assert(J >= 0 && J < 3);
  return base[I * 3 + J];
}

constexpr inline double& at(Mat3d& m, size_t i, size_t j) {
  double* base = &m.m11;
  return base[i * 3 + j];
}

template <typename T>
concept VecLike = requires(T v) {
  { v.x } -> std::convertible_to<double>;
  { v.y } -> std::convertible_to<double>;
  { v.z } -> std::convertible_to<double>;
};

} // namespace matrix_traits

/* ------------------------------------------------------------------------- */

class ParticleArray {
public:
  ParticleArray() {};

  inline void resize(size_t n) { raw.resize(3 * n); }
  inline double& operator()(size_t i, size_t j) { return raw.data()[i * 3 + j]; }
  inline double* data() { return raw.data(); }
  inline size_t size() { return raw.size(); }

private:
  std::vector<double> raw;
};

struct IOContext {

  struct IODataBuffer {
    size_t particle_count = 0;
    size_t type_count = 0;
    Mat3d cell{};
    Vec3d origin{};
    ParticleArray positions{};
    ParticleArray velocities{};
  };

  enum IOContextFlags : uint16_t {
    REMAP_ATOM  = 1 << 1,
    DOMAIN_ONLY = 1 << 2,
    TRICLINIC   = 1 << 3,
  };

  uint16_t flags = 0;
  IODataBuffer data{};

  template <uint16_t f>
  inline bool has() {
    return (flags & (f));
  }

  template<int16_t f>
  void set(bool toggle) {
    (toggle) ? flags |= f : flags &= f;
  }
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
      : m_path(std::move(path)),
        m_mode(mode),
        m_compression(compression) {}

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

template <typename T>
const Metadata get_parser_metadata() {
  return {};
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

  virtual inline bool operator()(IOContext&) { return false; };
  virtual inline bool operator()(IOContext&, size_t) { return false; };
  virtual size_t size() = 0;
};

class TextParser : public Parser {

protected:
  bool m_eof = false;
  std::vector<size_t> m_loc{};
  TextFile m_file;

  std::array<std::string_view, 10> m_line_buffer{};
  TokenSet m_tokens{};

  inline std::string_view& current_line(void) { return m_line_buffer[0]; };
  inline const std::string_view& current_line(void) const { return m_line_buffer[0]; }

public:
  TextParser(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression) {}

  virtual ~TextParser() override = default;

  virtual int64_t next() = 0;
  virtual bool parse(IOContext& ctx) = 0;

  bool operator()(IOContext& ctx) override;
  bool operator()(IOContext& ctx, size_t index) override;

  void scan();

  inline bool eof() { return m_eof; }
  inline size_t size() override { return m_loc.size(); };
};

/* ------------------------------------------------------------------------- */
using ParserCreator = std::function<std::unique_ptr<Parser>(std::string, FileMode, FileCompression)>;

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

  template <typename T>
  void register_format() {

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

template <bool RemapAtomID>
inline void assert_particle_index(size_t& index, const size_t id, const size_t n, const size_t particle_count) {
  if constexpr (RemapAtomID) {
    index = particle_count - 1;
  } else {
    index = cmp::le(id, n) ? (id - 1) : PANIC("ID=%ld is out-of-range N=%ld", id, n);
  }
}

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lmpdata {

using AtomDataParserFunction = bool (*)(const TokenSet&, IOContext&, size_t);
using KeywordParserFunction = bool (*)(const TokenSet&, IOContext&);

enum class AtomStyle {
  Atomic,
  Charge,
  Molecular,
  Full,
  Count,
  None,
};

constexpr TokenNeedles<enum_size<AtomStyle>()> NEEDLES_ATOM_STYLE{
    "atomic",
    "charge",
    "molecule",
    "full",
};

enum class Keyword {
  Atoms,
  AtomTypes,
  XAxis,
  YAxis,
  ZAxis,
  Tilts,
  Section,
  Count,
};

// Needles for header keywords (HDRKEY)
constexpr TokenNeedles<1> NEEDLES_KEYWORD_ATOMS{"atoms"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_ATOM_TYPES{"atom", "types"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_XAXIS{"xlo", "xhi"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_YAXIS{"ylo", "yhi"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_ZAXIS{"zlo", "zhi"};
constexpr TokenNeedles<3> NEEDLES_KEYWORD_TILTS{"xy", "xz", "yz"};

enum class Section {
  Masses,
  Atoms,
  Velocities,
  Bonds,
  Angles,
  Dihedrals,
  Count,
  Header,
  Ignored,
};

constexpr TokenNeedles<enum_size<Section>()> NEEDLES_SECTIONS{
    "Masses",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
};

bool parse_keyword_atoms(const TokenSet&, IOContext&);
bool parse_keyword_atom_types(const TokenSet&, IOContext&);
bool parse_keyword_xaxis(const TokenSet&, IOContext&);
bool parse_keyword_yaxis(const TokenSet&, IOContext&);
bool parse_keyword_zaxis(const TokenSet&, IOContext&);
bool parse_keyword_tilts(const TokenSet&, IOContext&);
bool parse_keyword_section(const TokenSet&, IOContext&);

constexpr std::array<KeywordParserFunction, enum_size<Keyword>()> keyword_parsers{
    &parse_keyword_atoms,
    &parse_keyword_atom_types,
    &parse_keyword_xaxis,
    &parse_keyword_yaxis,
    &parse_keyword_zaxis,
    &parse_keyword_tilts,
    &parse_keyword_section,
};


template <bool RemapAtomID>
inline bool parse_atom_data_atomic(const TokenSet& tokens, IOContext& ctx, size_t pcounter) {
  constexpr size_t len = 5; // atom-id atom-type x y z

  size_t index = 0;
  size_t id, type;
  double x, y, z;

  if (tokens.atleast(len) &&
      parse_tokens(
          tokens.at<0>(),
          id,
          tokens.at<1>(),
          type,
          tokens.at<2>(),
          x,
          tokens.at<3>(),
          y,
          tokens.at<4>(),
          z)) {

    assert_particle_index<RemapAtomID>(index, id, ctx.data.particle_count, pcounter);

    ctx.data.positions(index, 0) = x;
    ctx.data.positions(index, 1) = y;
    ctx.data.positions(index, 2) = z;

    ctx.data.velocities(index, 0) = 0.0;
    ctx.data.velocities(index, 1) = 0.0;
    ctx.data.velocities(index, 2) = 0.0;

    return true;
  }
  return false;
}

template <bool RemapAtomID>
inline bool parse_atom_data_charge(const TokenSet& tokens, IOContext& ctx, size_t pcounter) {
  constexpr size_t len = 6; // atom-id atom-type q x y z

  size_t index = 0;
  size_t id, type;
  double x, y, z;
  Parsable<double, false> q;

  if (tokens.atleast(len) &&
      parse_tokens(
          tokens.at<0>(),
          id,
          tokens.at<1>(),
          type,
          tokens.at<2>(),
          q,
          tokens.at<3>(),
          x,
          tokens.at<4>(),
          y,
          tokens.at<5>(),
          z)) {

    assert_particle_index<RemapAtomID>(index, id, ctx.data.particle_count, pcounter);

    ctx.data.positions(index, 0) = x;
    ctx.data.positions(index, 1) = y;
    ctx.data.positions(index, 2) = z;

    ctx.data.velocities(index, 0) = 0.0;
    ctx.data.velocities(index, 1) = 0.0;
    ctx.data.velocities(index, 2) = 0.0;

    return true;
  }
  return false;
}

template <bool RemapAtomID>
inline bool parse_atom_data_molecular(const TokenSet&, IOContext&, size_t) {
  [[maybe_unused]] constexpr size_t len = 6; // molecular: atom-id molecule-id atom-type x y z
  PANIC("Not Implemented style=molecular");
  return false;
}

template <bool RemapAtomID>
inline bool parse_atom_data_full(const TokenSet&, IOContext&, size_t) {
  [[maybe_unused]] constexpr size_t len = 7; // full:      atom-id molecule-id atom-type q x y z
  PANIC("Not Implemented style=full");
  return false;
}

template <bool RemapAtomID>
inline const std::array<const AtomDataParserFunction, enum_size<AtomStyle>()>& atom_data_parsers() {
  static constexpr std::array<const AtomDataParserFunction, enum_size<AtomStyle>()> fptr_hdlr{
      &parse_atom_data_atomic<RemapAtomID>,
      &parse_atom_data_charge<RemapAtomID>,
      &parse_atom_data_molecular<RemapAtomID>,
      &parse_atom_data_full<RemapAtomID>,
  };
  return fptr_hdlr;
}

template <bool RemapAtomID>
inline bool parse_atom_velocities(const TokenSet& tokens, IOContext& ctx, size_t pcounter) {
  constexpr size_t len = 4;

  size_t id, index = 0;
  double vx, vy, vz;

  if (tokens.atleast(len) &&
      parse_tokens(tokens.at<0>(), id, tokens.at<1>(), vx, tokens.at<2>(), vy, tokens.at<3>(), vz)) {

    assert_particle_index<RemapAtomID>(index, id, ctx.data.particle_count, pcounter);

    ctx.data.velocities(index, 0) = vx;
    ctx.data.velocities(index, 1) = vy;
    ctx.data.velocities(index, 2) = vz;

    return true;
  }
  return false;
}

class LAMMPSDataParser : public TextParser {
public:
  LAMMPSDataParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression),
        m_atom_style(AtomStyle::None),
        m_current_section(Section::Header) {}

  int64_t next() override;
  bool parse(IOContext& proxy) override;

private:
  AtomStyle m_atom_style;
  Section m_current_section;

  // keep track of the current section line
  inline std::string_view& section_line(void) { return m_line_buffer[1]; }
  inline const std::string_view& section_line(void) const { return m_line_buffer[1]; }

  bool parse_section_header(IOContext& ctx);
  bool parse_section_atoms(IOContext& ctx);
  bool parse_section_masses(void);
  bool parse_section_velocities(IOContext& ctx);
  void jumpto_next_section();

  inline void set_current_section(std::string_view str) {
    size_t index;
    if (token_match_any(str, NEEDLES_SECTIONS, index)) {
      m_current_section = enum_from_index<Section>(index);
    } else {
      m_current_section = Section::Ignored;
    }
  }
};

} // namespace lmpdata

template <>
inline const Metadata get_parser_metadata<lmpdata::LAMMPSDataParser>() {
  return {
      .name = "LAMMPS Data",
      .aliases = {"lmp", "lmp-data", "data"},
      .description = "LAMMPS Data format",
  };
}

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Dump format

namespace lmpdump {

enum class Field {
  ID,      // atom id
  Type,    // atom type
  Element, // element
  X,       // x wrapped
  Y,       // y wrapped
  Z,       // z wrapped
  XU,      // x unwrapped
  YU,      // y unwrapped
  ZU,      // z unwrapped
  XS,      // x scaled position
  YS,      // y scaled position
  ZS,      // z scaled position
  XSU,     // x scaled unwrapped position
  YSU,     // y scaled unwrapped position
  ZSU,     // z scaled unwrapped position
  VX,      // x velocity
  VY,      // y velocity
  VZ,      // z velocity
  None,    // ignored field
  Count,   //
};

constexpr TokenNeedles<enum_size<Field>()> NEEDLES_FIELDS = {
    "id",
    "type",
    "element",
    "x",
    "y",
    "z",
    "xu",
    "yu",
    "zu",
    "xs",
    "ys",
    "zs",
    "xsu",
    "ysu",
    "zsu",
    "vx",
    "vy",
    "vz",
};

struct Coordinates {
  bool scaled;
  bool unwrapped;
  Field x, y, z;
};

constexpr std::array<Coordinates, 4> COORDINATES{{
    {false, false, Field::X, Field::Y, Field::Z},
    {false, true, Field::XU, Field::YU, Field::ZU},
    {true, false, Field::XS, Field::YS, Field::ZS},
    {true, true, Field::XSU, Field::YS, Field::ZS},
}};

struct Property {
  Field field = Field::None;
  size_t index = 0;
};


struct Properties : properties_traits::GenericContainer<Property, Field, &Property::field> {
  bool triclinic = false;
  bool scaled = false;
  bool unwrapped = false;
  IJK pos_ijk{}, vel_ijk{};

  bool (*get_pos)(Vec3d&, const TokenSet&, const IJK&, const IOContext&) = nullptr;
  bool (*get_vel)(Vec3d&, const TokenSet&, const IJK&) = nullptr;
  bool (*get_id)(size_t&, const TokenSet&, size_t, size_t) = nullptr;
  bool (*get_type)(size_t&, const TokenSet&, size_t) = nullptr;

  template <Field i, Field j, Field k>
  inline constexpr IJK ijk() const {
    return IJK{get<i>().index, get<j>().index, get<k>().index};
  }

  inline IJK ijk(Field i, Field j, Field k) { return {get(i).index, get(j).index, get(k).index}; }

  inline bool has_positions() {
    return (has<Field::X, Field::Y, Field::Z>() || has<Field::XU, Field::YU, Field::ZU>() ||
            has<Field::XS, Field::YS, Field::ZS>() || has<Field::XSU, Field::YSU, Field::ZSU>());
  }
};

template <bool IsTriclinic, matrix_traits::VecLike T>
inline bool parse_box_bounds_axis(const TokenSet& tokens, const T& v) {
  if constexpr (IsTriclinic) {
    return tokens.atleast(3) && parse_tokens(tokens.at<0>(), v.x, tokens.at<1>(), v.y, tokens.at<2>(), v.z);
  } else {
    return tokens.atleast(2) && parse_tokens(tokens.at<0>(), v.x, tokens.at<1>(), v.y);
  }
}

inline void convert_bbox_to_lattice_vectors(IOContext& ctx, const Mat3d& bbox) {
  double xlo = bbox.m11 - std::min(bbox.m13 + bbox.m23, std::min(0., bbox.m23));
  double xhi = bbox.m12 - std::max(bbox.m13 + bbox.m23, std::max(0., bbox.m23));
  double ylo = bbox.m21 - std::min(0., bbox.m33);
  double yhi = bbox.m22 - std::max(0., bbox.m33);
  double zlo = bbox.m31;
  double zhi = bbox.m32;

  ctx.data.cell.m11 = xhi - xlo;
  ctx.data.cell.m22 = yhi - ylo;
  ctx.data.cell.m33 = zhi - zlo;
  ctx.data.cell.m12 = bbox.m13;
  ctx.data.cell.m13 = bbox.m23;
  ctx.data.cell.m23 = bbox.m33;

  ctx.data.origin.x = xlo;
  ctx.data.origin.y = ylo;
  ctx.data.origin.z = zlo;
};

template <bool HasVelocities>
inline bool parse_velocities(Vec3d& v, const TokenSet& tokens, const IJK& inds) {
  if constexpr (HasVelocities) {
    return parse_tokens(tokens[inds.i], v.x, tokens[inds.j], v.y, tokens[inds.k], v.z);
  } else {
    v = {0., 0., 0};
    return true;
  }
};

template <bool IsScaled, bool IsUnwrapped, bool IsTriclinic>
inline bool parse_positions(Vec3d& r, const TokenSet& tokens, const IJK& inds, const IOContext& ctx) {
  if (!parse_tokens(tokens[inds.i], r.x, tokens[inds.j], r.y, tokens[inds.k], r.z)) {
    return false;
  }

  if constexpr (IsScaled) {
    if constexpr (IsTriclinic) {
      // x = xlo + x * lx + y * xy + z * xz
      // y = zlo + y * ly + z * xz
      // z = ylo + z * lz;
      r.x = ctx.data.origin.x + r.x * ctx.data.cell.m11 + r.y + ctx.data.cell.m12 + r.z * ctx.data.cell.m13;
      r.y = ctx.data.origin.y + r.y * ctx.data.cell.m22 + r.z * ctx.data.cell.m13;
      r.z = ctx.data.origin.z + r.z * ctx.data.cell.m33;
    } else {
      r.x = ctx.data.origin.x + r.x * ctx.data.cell.m11;
      r.y = ctx.data.origin.y + r.y * ctx.data.cell.m22;
      r.z = ctx.data.origin.z + r.z * ctx.data.cell.m33;
    }
  }

  return true;
}

template <bool HasID>
inline bool parse_atom_id(size_t& id, const TokenSet& tokens, size_t index, size_t particle_count) {
  if constexpr (HasID) {
    return parse_tokens(tokens[index], id);
  } else {
    id = particle_count;
    return true;
  }
}

template <bool HasParticleType>
inline bool parse_atom_type(size_t& type, const TokenSet& tokens, size_t index) {
  if constexpr (HasParticleType) {
    return parse_tokens(tokens[index], type);
  } else {
    type = 0; // set default type to 0
    return true;
  }
}

template <bool IsTriclinic>
inline bool parse_atom_properties(const TokenSet& tokens, Properties& properties) {

  // parse ITEM: ATOMS ....
  size_t ifield, headlen = 2;
  for (size_t i = headlen; i < tokens.size(); ++i) {
    if (token_match_any(tokens[i], NEEDLES_FIELDS, ifield)) {
      properties.set({.field = enum_from_index<Field>(ifield), .index = i - headlen});
    }
  }

  properties.get_id = (properties.has<Field::ID>()) ? parse_atom_id<true> : parse_atom_id<false>;
  properties.get_type = (properties.has<Field::Type>()) ? parse_atom_type<true> : parse_atom_type<false>;

  // deduced coodinate system base on present field
  for (const Coordinates& c : COORDINATES) {
    if (properties.has(c.x, c.y, c.z)) {
      properties.scaled = c.scaled;
      properties.unwrapped = c.unwrapped;
      properties.pos_ijk = properties.ijk(c.x, c.y, c.z);
      properties.expect(c.x, c.y, c.z);
    }
  }

  if (!properties.has_positions()) {
    PANIC("No position fields are defined");
  }

  properties.get_pos =
      (properties.scaled)
          ? parse_positions<true, false, IsTriclinic>
          : parse_positions<false, false, IsTriclinic>;

  // check if there are velocity fields
  if (properties.has<Field::VX, Field::VY, Field::VZ>()) {
    properties.get_vel = parse_velocities<true>;
    properties.vel_ijk = properties.ijk<Field::VX, Field::VY, Field::VZ>();
    properties.expect(Field::VX, Field::VY, Field::VZ);
  } else {
    properties.get_vel = parse_velocities<false>;
  }

  return true;
}

template <bool RemapAtomID>
inline bool parse_atom_data(const TokenSet& tokens, IOContext& ctx, const Properties& p, size_t particle_count) {

  size_t id, type;
  Vec3d pos, vel;
  
  if (!(tokens.atleast(p.len()) && p.get_id(id, tokens, p.get<Field::ID>().index, particle_count) &&
        p.get_pos(pos, tokens, p.pos_ijk, ctx) && p.get_vel(vel, tokens, p.vel_ijk) &&
        p.get_type(type, tokens, p.get<Field::Type>().index))) {
    return false;
  }

  size_t particle_index;
  assert_particle_index<RemapAtomID>(particle_index, id, ctx.data.particle_count, particle_count);

  ctx.data.positions(particle_index, 0) = pos.x;
  ctx.data.positions(particle_index, 1) = pos.y;
  ctx.data.positions(particle_index, 2) = pos.z;
  
  ctx.data.velocities(particle_index, 0) = vel.x;
  ctx.data.velocities(particle_index, 1) = vel.y;
  ctx.data.velocities(particle_index, 2) = vel.z;

  return true;
}

class LAMMPSDumpParser : public TextParser {
public:
  LAMMPSDumpParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;
};

} // namespace lmpdump

template <>
inline const Metadata get_parser_metadata<lmpdump::LAMMPSDumpParser>() {
  return {
      .name = "LAMMPS Dump",
      .aliases = {"dump", "lmp-dump", "lammps-dump"},
      .description = "LAMMPS Dump format",
  };
}

/* ------------------------------------------------------------------------- */
// Parser for Extended XYZ Data format

namespace xyz {

enum class FieldType : uint8_t {
  R = 1 << 0,
  S = 1 << 1,
  I = 1 << 2,
  Count = 3,
  None = 0,
};

constexpr TokenNeedles<enum_size<FieldType>()> NEEDLES_FIELD_TYPE{
    "R",
    "S",
    "I",
};

enum class Field {
  Species,
  Positions,
  Velocities,
  ID,
  Type,
  None,
  Count,
};

constexpr TokenNeedles<enum_size<Field>()> NEEDLES_FIELDS{
    "species",
    "pos",
    "velo",
    "id",
    "type",
};

struct Property {
  Field field = Field::None;
  FieldType type = FieldType::None;
  size_t icol = 0;
  size_t ncol = 0;
};

struct Properties : public properties_traits::GenericContainer<Property, Field, &Property::field> {

  using Enum = typename properties_traits::GenericContainer<Property, Field, &Property::field>::enum_t;

  template <Enum... args>
  constexpr void expect() {
    // static_assert((std::is_same_v<Enum, Enum> && args...), "Only E enum type allowed");
    ((length += (has<args>()) ? get<args>().ncol : static_cast<size_t>(0)), ...);
  }

  bool (*get_id)(size_t&, const TokenSet&, size_t, size_t) = nullptr;
  bool (*get_type)(size_t&, const TokenSet&, size_t) = nullptr;
  bool (*get_vel)(Vec3d&, const TokenSet&, const IJK&) = nullptr;

  IJK pos_ijk{}, vel_ijk{};
};

template <size_t N, size_t M>
bool parse_xyz_comment_line(StringView str, TokenSetTmpl<N>& lattice, TokenSetTmpl<M>& properties) {
  constexpr StringView key_lattice = "Lattice=\"";
  constexpr size_t key_lattice_len = key_lattice.size(); // 9
  constexpr StringView key_properties = "Properties=";
  constexpr size_t key_properties_len = key_properties.size(); // 11

  const char* begin = str.data();
  const char* end = str.data() + str.size();

  if (!startswith(str, key_lattice)) {
    return false;
  }

  // 1. Parse Lattice=...
  const char* first = begin + key_lattice_len;
  const char* last = first;

  while (last < end && !DoubleQuoteDelimiter{}(last)) {
    ++last;
  }

  if (last >= end) {
    return false;
  }

  tokenize(StringView(first, last), lattice);

  // 2. Parse Properties=....
  first = last;
  while (first + key_properties_len <= end) {

    if (std::memcmp(first, key_properties.data(), key_properties_len) == 0) {
      first += key_properties_len;
      last = first;

      while (last < end && !SpaceDelimiter{}(last)) {
        ++last;
      }

      tokenize_impl(StringView(first, last), properties, ColumnDelimiter{});
    }

    ++first;
  }

  return true;
}

inline bool parse_xyz_lattice(const TokenSet& tokens, IOContext& ctx) {
  return tokens.atleast(9) &&
         parse_tokens(
             tokens.at<0>(),
             ctx.data.cell.m11, // ax
             tokens.at<1>(),
             ctx.data.cell.m21, // ay
             tokens.at<2>(),
             ctx.data.cell.m31, // az
             tokens.at<3>(),
             ctx.data.cell.m12, // bx
             tokens.at<4>(),
             ctx.data.cell.m22, // by
             tokens.at<5>(),
             ctx.data.cell.m32, // bz
             tokens.at<6>(),
             ctx.data.cell.m13, // cx
             tokens.at<7>(),
             ctx.data.cell.m23, // cy
             tokens.at<8>(),
             ctx.data.cell.m33  // cz
         );
}

template <bool HasID>
inline bool parse_atom_id(size_t& id, const TokenSet& tokens, size_t index, size_t particle_count) {
  if constexpr (HasID) {
    return parse_tokens(tokens[index], id);
  } else {
    id = particle_count;
    return true;
  }
}

template <bool HasParticleType>
inline bool parse_atom_type(size_t& type, const TokenSet& tokens, size_t index) {
  if constexpr (HasParticleType) {
    return parse_tokens(tokens[index], type);
  } else {
    type = 0; // set default type to 0
    return true;
  }
}

template <bool HasVelocities>
inline bool parse_velocities(Vec3d& v, const TokenSet& tokens, const IJK& inds) {
  if constexpr (HasVelocities) {
    return parse_tokens(tokens[inds.i], v.x, tokens[inds.j], v.y, tokens[inds.k], v.z);
  } else {
    v = {0., 0., 0};
    return true;
  }
}

inline bool parse_positions(Vec3d& r, const TokenSet& tokens, const IJK& inds) {
  return parse_tokens(tokens[inds.i], r.x, tokens[inds.j], r.y, tokens[inds.k], r.z);
}

inline bool parse_xyz_properties(const TokenSet& tokens, Properties& properties) {

  if (tokens.size() % 3 != 0) {
    return false;
  }

  size_t icol = 0, ncol, ifield, itype;

  properties.clear();

  for (size_t i = 0; i < tokens.size(); i += 3) {
    size_t j = i + 1;
    size_t k = i + 2;

    parse_tokens(tokens[k], ncol);

    if ((token_match_any(tokens[i], NEEDLES_FIELDS, ifield) &&
         token_match_any(tokens[j], NEEDLES_FIELD_TYPE, itype))) {
      properties.set({
          .field = enum_from_index<Field>(ifield),
          .type = enum_from_index<FieldType>(itype),
          .icol = icol,
          .ncol = ncol,
      });
    }

    icol += ncol;
  }

  // Field validation
  if (!properties.has<Field::Positions>()) {
    PANIC("Positions are not defined in properties");
  }
  size_t ipos = properties.get<Field::Positions>().icol;
  properties.pos_ijk = {ipos, ipos + 1, ipos + 2};

  if (!(properties.has<Field::Type>() || properties.has<Field::Species>())) {
    PANIC("Either Species or TYPE should be define in xyz properties")
  }

  properties.expect<Field::Species, Field::Positions, Field::Velocities, Field::ID, Field::Type>();
  properties.get_id = (properties.has<Field::ID>()) ? parse_atom_id<true> : parse_atom_id<false>;
  properties.get_type = (properties.has<Field::Type>()) ? parse_atom_type<true> : parse_atom_type<false>;

  if (properties.has<Field::Velocities>()) {
    properties.get_vel = parse_velocities<true>;
    size_t ivel = properties.get<Field::Velocities>().icol;
    properties.vel_ijk = {ivel, ivel + 1, ivel + 2};
  } else {
    properties.get_vel = parse_velocities<false>;
  }

  return true;
}

template<bool RemapAtomID>
inline bool parse_xyz_atom_line(
    const TokenSet& tokens,
    IOContext& ctx,
    const Properties& properties,
    size_t particle_count) {

  size_t id, particle_index;
  Vec3d r, v;

  if (!(tokens.atleast(properties.len()) && parse_positions(r, tokens, properties.pos_ijk) &&
        properties.get_vel(v, tokens, properties.vel_ijk))) {
    return false;
  }

  if constexpr (RemapAtomID) {
    properties.get_id(id, tokens, properties.get<Field::ID>().icol, particle_count);
  } else {
    id = particle_count;
  }

  assert_particle_index<RemapAtomID>(particle_index, id, ctx.data.particle_count, particle_count);

  ctx.data.positions(particle_index, 0) = r.x;
  ctx.data.positions(particle_index, 1) = r.y;
  ctx.data.positions(particle_index, 2) = r.z;

  ctx.data.velocities(particle_index, 0) = v.x;
  ctx.data.velocities(particle_index, 1) = v.y;
  ctx.data.velocities(particle_index, 2) = v.z;

  return true;
};

class ExtendedXYZParser : public TextParser {
public:
  ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;

private:
  TokenSet m_ptokens;
  Properties m_properties;
};

} // namespace xyz

// clang-format off
template <> inline const Metadata get_parser_metadata<xyz::ExtendedXYZParser>() {
  return {
    .name        = "Extended XYZ",
    .aliases     = { "xyz", "exyz", "ext-xyz" },
    .description = "Extended XYZ"
  };
}
// clang-format on

/* ------------------------------------------------------------------------- */

struct FileInfos {
  FileCompression compression = FileCompression::NONE;
  const RegisteredFormat format{};
};

std::string remove_leading_chars(const std::string& str, const char c);
const FileInfos guess_file_infos(const std::string& filepath);
const FileInfos get_file_infos(const std::string& filepath, std::string format, std::string compression);

inline void read_atom_file(
    IOContext& ctx,
    const std::string& file,
    const std::string& format = "",
    const std::string& compression = "") {

  // Ensure the file exists
  std::string filepath = std::string(std::filesystem::absolute(file));
  if (!std::filesystem::exists(filepath)) {
    PANIC("Input file doest not exists: '%s'\n", filepath.c_str());
  }

  // create reader object
  FileInfos infos = get_file_infos(filepath, format, compression);
  std::unique_ptr<Parser> parser = infos.format.creator(filepath, FileMode::READ, infos.compression);
  if (parser == nullptr) {
    PANIC("Something went wrong... I AM PANICKING");
  }

  if (!(*parser)(ctx)) {
    PANIC("Something went wrong when parsing the file");
  }
}

} // namespace io
