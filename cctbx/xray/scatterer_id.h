// (c) O.V.D., OlexSys Ltd, 2025
#pragma once
#include <limits>
#include <charconv>
#include <bit>
#include <array>
#include <iomanip>

namespace cctbx {
  namespace xray {
    /* Stores coordinates with float type precision, adopted from Olex2. */

    template <typename FloatType, class crd_t, class mask_info, uint64_t cell_m>
    struct scatterer_id_base : public mask_info {
      uint64_t id;

      scatterer_id_base(const scatterer_id_base& sid)
        : id(sid.id)
      {}

      scatterer_id_base(uint64_t id = ~0)
        : id(id)
      {}

      scatterer_id_base(std::istream& is)
      {
        is.read((char*)&id, sizeof(uint64_t));
      }

      scatterer_id_base(std::string id_string)
      {
        std::stringstream ss;
        ss << std::hex << id_string;
        ss >> id;
      }

      // use multiplier like 1e3 to 'crop' values
      scatterer_id_base(int8_t z, const crd_t& crd, short data,
        FloatType multiplier = 1)
      {
        id = ((uint64_t)z) & mask_info::z_mask;
        static const int64_t k = mask_info::mask_m / cell_m;
        int64_t x = multiplier == 1 ? (int64_t)(crd[0] * k)
          : ((int64_t)round(crd[0] * multiplier)) / multiplier * k;
        if (x < 0) {
          id |= mask_info::a_sig;
          id |= (((-x) << 8) & mask_info::a_mask);
        }
        else {
          id |= ((std::abs(x) << 8) & mask_info::a_mask);
        }
        x = multiplier == 1 ? (int64_t)(crd[1] * k)
          : ((int64_t)round(crd[1] * multiplier)) / multiplier * k;
        if (x < 0) {
          id |= mask_info::b_sig;
          id |= (((-x) << (8 + mask_info::a_shift)) & mask_info::b_mask);
        }
        else {
          id |= ((x << (8 + mask_info::a_shift)) & mask_info::b_mask);
        }
        x = multiplier == 1 ? (int64_t)(crd[2] * k)
          : ((int64_t)round(crd[2] * multiplier)) / multiplier * k;
        if (x < 0) {
          id |= mask_info::c_sig;
          id |= (((-x) << (8 + mask_info::a_shift * 2)) & mask_info::c_mask);
        }
        else {
          id |= ((x << (8 + mask_info::a_shift * 2)) & mask_info::c_mask);
        }
        x = data;
        id |= (((int64_t)data << (8 + mask_info::a_shift * 3) + 3) & mask_info::d_mask);
      }

      bool test(const crd_t& crd) const {
        return std::abs(crd[0]) <= cell_m &&
          std::abs(crd[1]) <= cell_m &&
          std::abs(crd[2]) <= cell_m;
      }

      scatterer_id_base& operator = (const scatterer_id_base& i) {
        this->id = i.id;
        return *this;
      }

      scatterer_id_base& operator = (const uint64_t& i) {
        this->id = i;
        return *this;
      }

      bool operator == (const scatterer_id_base& i) const {
        return id == i.id;
      }

      bool operator != (const scatterer_id_base& i) const {
        return id != i.id;
      }

      int Compare(const scatterer_id_base& i) const {
        return id < i.id ? -1 : (id > i.id ? 1 : 0);
      }

      crd_t get_crd() const {
        static const FloatType k = (FloatType)cell_m / mask_info::mask_m;
        crd_t r(
          static_cast<FloatType>((id & mask_info::a_mask) >> 8) * k,
          static_cast<FloatType>((id & mask_info::b_mask) >> (8 + mask_info::a_shift)) * k,
          static_cast<FloatType>((id & mask_info::c_mask) >> (8 + mask_info::a_shift * 2)) * k);
        if ((id & mask_info::a_sig) != 0) {
          r[0] = -r[0];
        }
        if ((id & mask_info::b_sig) != 0) {
          r[1] = -r[1];
        }
        if ((id & mask_info::c_sig) != 0) {
          r[2] = -r[2];
        }
        return r;
      }

      int8_t get_z() const {
        return (int8_t)(id & mask_info::z_mask);
      }

      // this wil always be returned as unsigned value
      short get_data() const {
        return (short)((id & mask_info::d_mask) >> (8 + mask_info::a_shift * 3 + 3));
      }

    };

    /*  0-8 - z, 8-25, 25-42, 42-59 - a, b c, 59-61 - signs, 62-64 - data
    precision 16: ~0.0025, 1: ~0.0000077
    */
    struct scatterer_id_masks_d2 {
      const static uint64_t
        z_mask = 0x00000000000000FF,
        a_mask = 0x0000000001FFFF00,
        b_mask = 0x000003FFFE000000,
        c_mask = 0x07FFFC0000000000,
        a_sig  = 0x0800000000000000,
        b_sig  = 0x1000000000000000,
        c_sig  = 0x2000000000000000,
        d_mask = 0xC000000000000000,
        mask_m = 0x000000000001FFFF; // max crd value
      const static int a_shift = 17;
    };

    /* 2 extra bits for data */
    template <typename FloatType, class crd_t, uint64_t cell_m>
    struct scatterer_id_2
      : public scatterer_id_base<FloatType, crd_t, scatterer_id_masks_d2, cell_m> {
      typedef scatterer_id_base<FloatType, crd_t, scatterer_id_masks_d2, cell_m> parent_t;
      scatterer_id_2(const scatterer_id_2& sid)
        : parent_t(sid)
      {}

      scatterer_id_2(uint64_t id = ~0)
        : parent_t(id)
      {}

      scatterer_id_2(int8_t z, const crd_t& crd, short data,
        FloatType multiplier = 1)
        : parent_t(z, crd, data, multiplier)
      {}

      scatterer_id_2(std::string id_string)
       :parent_t(id_string)
      { }
    };

    /*  0-8 - z, 8-24, 24-40, 40-56 - a, b c, 56-59 - signs, 59-64 - data
    precision 16: ~0.004, 1: ~0.00002
    */
    struct scatterer_id_masks_d5 {
      const static uint64_t
        z_mask = 0x00000000000000FF,
        a_mask = 0x0000000000FFFF00,
        b_mask = 0x000000FFFF000000,
        c_mask = 0x00FFFF0000000000,
        a_sig =  0x0100000000000000,
        b_sig =  0x0200000000000000,
        c_sig =  0x0400000000000000,
        d_mask = 0xF800000000000000,
        mask_m = 0x000000000000FFFF; // max crd value
      const static int a_shift = 16;
    };

    /* 5 extra bits for data */
    template <typename FloatType, class crd_t, uint64_t cell_m>
    struct scatterer_id_5
      : public scatterer_id_base<FloatType, crd_t, scatterer_id_masks_d5, cell_m> {
      typedef scatterer_id_base<FloatType, crd_t, scatterer_id_masks_d5, cell_m> parent_t;
      scatterer_id_5(const scatterer_id_5& sid)
        : parent_t(sid)
      {}

      scatterer_id_5(uint64_t id = ~0)
        : parent_t(id)
      {}
      // use multiplier like 1e3 to 'crop' values
      scatterer_id_5(int8_t z, const crd_t& crd, short data,
        FloatType multiplier = 1)
        : parent_t(z, crd, data, multiplier)
      {}

      scatterer_id_5(std::string id_string)
       : parent_t(id_string)
      { }
    };

    namespace atomID_utils {

        inline constexpr double max_fractional_coordinate = 16.0;

        inline const double integer_max = static_cast<double>(std::numeric_limits<int32_t>::max());

        inline const double encode_scale = integer_max / max_fractional_coordinate;

        inline const double decode_scale = max_fractional_coordinate / integer_max;

    } // namespace atomID_utils
    template <typename FloatType, class crd_t>
    class scatterer_id_big {
    private:
        int32_t frac_x_int;
        int32_t frac_y_int;
        int32_t frac_z_int;
        int16_t data_;
        uint8_t Z_;
        uint8_t reserved_;

        static std::int32_t encode_coordinate(FloatType coordinate)
        {
            if (!std::isfinite(coordinate)) {
                throw std::invalid_argument("Coordinate must be finite");
            }

            if (coordinate < -atomID_utils::max_fractional_coordinate ||
                coordinate > atomID_utils::max_fractional_coordinate) {
                throw std::out_of_range(
                    "Coordinate must be between -16 and 16"
                );
            }

            return static_cast<std::int32_t>(
                std::llround(coordinate * atomID_utils::encode_scale)
            );
        }

        static FloatType decode_coordinate(std::int32_t encoded)
        {
            return static_cast<FloatType>(encoded) * atomID_utils::decode_scale;
        }

        void validate_loaded_data() const
        {
            if (Z_ == 0) {
                throw std::runtime_error(
                    "Invalid atom ID: atomic number is zero");
            }

            /*
            * The encoder maps -16 to -INT32_MAX, not INT32_MIN.
            * Therefore, INT32_MIN cannot be produced by a valid encoder.
            */
            if (frac_x_int == std::numeric_limits<std::int32_t>::min() ||
                frac_y_int == std::numeric_limits<std::int32_t>::min() ||
                frac_z_int == std::numeric_limits<std::int32_t>::min()) {
                throw std::runtime_error(
                    "Invalid atom ID: encoded coordinate is out of range");
            }

            /*
            * If reserved_ is currently required to be zero, enable this:
            *
            * if (reserved_ != 0) {
            *     throw std::runtime_error(
            *         "Unsupported atom ID format version");
            * }
            */
        }

        static scatterer_id_big<FloatType, crd_t> from_uint64(
            const std::uint64_t first,
            const std::uint64_t second) noexcept
        {
            scatterer_id_big<FloatType, crd_t> result;

            const auto frac_x_bits =
                static_cast<std::uint32_t>(first);

            const auto frac_y_bits =
                static_cast<std::uint32_t>(first >> 32);

            const auto frac_z_bits =
                static_cast<std::uint32_t>(second);

            const auto data_bits =
                static_cast<std::uint16_t>(second >> 32);

            result.frac_x_int = std::bit_cast<std::int32_t>(frac_x_bits);
            result.frac_y_int = std::bit_cast<std::int32_t>(frac_y_bits);
            result.frac_z_int = std::bit_cast<std::int32_t>(frac_z_bits);
            result.data_ = std::bit_cast<std::int16_t>(data_bits);

            result.Z_ =
                static_cast<std::uint8_t>(second >> 48);

            result.reserved_ =
                static_cast<std::uint8_t>(second >> 56);

            return result;
        }

        static scatterer_id_big<FloatType, crd_t> from_hex_string(std::string hex)
        {
            if (hex.size() != 32) {
                throw std::invalid_argument(
                    "scatterer_id_big hexadecimal string must contain exactly 32 characters");
            }

            std::uint64_t first{};
            std::uint64_t second{};

            const auto first_part = hex.substr(0, 16);
            const auto second_part = hex.substr(16, 16);

            const auto first_result = std::from_chars(
                first_part.data(),
                first_part.data() + first_part.size(),
                first,
                16);

            const auto second_result = std::from_chars(
                second_part.data(),
                second_part.data() + second_part.size(),
                second,
                16);

            if (first_result.ec != std::errc{} ||
                first_result.ptr != first_part.data() + first_part.size() ||
                second_result.ec != std::errc{} ||
                second_result.ptr != second_part.data() + second_part.size()) {
                throw std::invalid_argument(
                    "scatterer_id_big string contains invalid hexadecimal characters");
            }
            return from_uint64(first, second);
        }
    
    public:
        scatterer_id_big() : frac_x_int(0), frac_y_int(0), frac_z_int(0), data_(0), Z_(0), reserved_(0) {}

        scatterer_id_big(const crd_t& site, const int data, const int Z, const int reserved = 0){
            *this = scatterer_id_big(site[0], site[1], site[2], data, Z, reserved);
        }

        scatterer_id_big(const FloatType frac_x, const FloatType frac_y, const FloatType frac_z, const int data, const int Z, const int reserved = 0)
        {
            if (Z < 1 || Z > 255) {
                throw std::out_of_range("Z must be between 1 and 255");
            }
            if (data < std::numeric_limits<std::int16_t>::min() ||
                data > std::numeric_limits<std::int16_t>::max()) {
                throw std::out_of_range("data does not fit into int16_t");
            }
            if (reserved < 0 || reserved > 255) {
                throw std::out_of_range(
                    "reserved must be between 0 and 255");
            }

            frac_x_int = encode_coordinate(frac_x);
            frac_y_int = encode_coordinate(frac_y);
            frac_z_int = encode_coordinate(frac_z);
            data_ = static_cast<int16_t>(data);
            Z_ = static_cast<uint8_t>(Z);
            reserved_ = static_cast<uint8_t>(reserved);
        }

        scatterer_id_big(std::istream& input) {
            input.read(reinterpret_cast<char*>(this), sizeof(scatterer_id_big));
            if (!input) {
                throw std::runtime_error("Failed to read scatterer ID from stream");
            }
            validate_loaded_data();
        }

        scatterer_id_big(std::string id_string)
        {
            *this = from_hex_string(id_string);
            validate_loaded_data();
        }

        std::pair<uint64_t,uint64_t> as_uint64() const noexcept
        {
            // First 64-bit integer:
            // bits  0-31: frac_x_int
            // bits 32-63: frac_y_int
            const std::uint64_t first =
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(frac_x_int))
                |
                (static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(frac_y_int)) << 32);

            // Second 64-bit integer:
            // bits  0-31: frac_z_int
            // bits 32-47: data_
            // bits 48-55: Z_
            // bits 56-63: reserved_
            const std::uint64_t second =
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(frac_z_int))
                |
                (static_cast<std::uint64_t>(
                    static_cast<std::uint16_t>(data_)) << 32)
                |
                (static_cast<std::uint64_t>(Z_) << 48)
                |
                (static_cast<std::uint64_t>(reserved_) << 56);

            return { first, second };
        }

        std::string to_hex_string() const
        {
            const auto [first, second] = as_uint64();

            std::ostringstream stream;
            stream << std::hex
                << std::setfill('0')
                << std::setw(16) << first
                << std::setw(16) << second;

            return stream.str();
        }
        // Equality operator
        bool operator==(const scatterer_id_big& other) const {
            return frac_x_int == other.frac_x_int &&
                frac_y_int == other.frac_y_int &&
                frac_z_int == other.frac_z_int &&
                data_ == other.data_ &&
                Z_ == other.Z_ &&
                reserved_ == other.reserved_;
        }

        const bool is_initialized() const {
            return Z_ != 0;
        }

        std::string to_bytes() const {
            if (!is_initialized()) {
                throw std::runtime_error("Atom ID is not initialized");
            }
            std::ostringstream oss;
            oss.write(reinterpret_cast<const char*>(this), sizeof(scatterer_id_big));
            return oss.str();
        }

        crd_t get_crd() const {
            crd_t r(
                decode_coordinate(frac_x_int),
                decode_coordinate(frac_y_int),
                decode_coordinate(frac_z_int)
            );
            return r;
        }

        uint8_t get_z() const {
            return Z_;
        }

        std::int16_t get_data() const {
            return data_;
        }

        std::uint8_t reserved() const noexcept
        {
            return reserved_;
        }
    };

}} // namespace cctbx::xray
