#include "share/field/field_identifier.hpp"
#include "ekat/util/ekat_string_utils.hpp"

namespace scream
{

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const Units& units,
                 const std::string& grid_name)
 : FieldIdentifier(name,layout,units,grid_name,"real")
{
  // Nothing to do here
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const Units& units,
                 const std::string& grid_name,
                 const std::string& data_type)
 : m_name      (name)
 , m_units     (units)
 , m_grid_name (grid_name)
{
  // The list of supported data type does not contain 'Real'.
  // If data_type=="real", we find out what Real maps to, and use that.
  // That's because ekat::TypeMap does not allow duplicates in the keys,
  // and Real is just an alias to either float or double.
  m_data_type = data_type=="real" ? field_valid_data_types().at<Real>() : data_type;

  EKAT_REQUIRE_MSG (is_valid_field_data_type(m_data_type),
      "Error! Unsupported data type '" + data_type + "'.\n");
  set_layout (layout);
}

void FieldIdentifier::set_layout (const layout_type& layout) {
  set_layout(std::make_shared<layout_type>(layout));
}

void FieldIdentifier::set_layout (const layout_ptr_type& layout) {
  EKAT_REQUIRE_MSG (!m_layout,
      "Error! You cannot reset the layout once it's set.\n");
  EKAT_REQUIRE_MSG (layout,
      "Error! Invalid input layout pointer.\n");
  EKAT_REQUIRE_MSG (layout->are_dimensions_set(),
      "Error! Input layout must have dimensions set.\n");

  m_layout = layout;
  update_identifier ();
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "[" + m_grid_name + "] <" + m_data_type;
  if (m_layout->rank()>0) {
    m_identifier += ":" + e2str(m_layout->tags()[0]);
    for (int dim=1; dim<m_layout->rank(); ++dim) {
      m_identifier += "," + e2str(m_layout->tags()[dim]);
    }
    m_identifier += ">(" + std::to_string(m_layout->dims()[0]);
    for (int dim=1; dim<m_layout->rank(); ++dim) {
      m_identifier += "," + std::to_string(m_layout->dims()[dim]);
    }
    m_identifier += ")";
  } else {
    m_identifier += ">";
  }

  m_identifier += " [" + m_units.get_string() + "]";
}

// Free functions for identifiers comparison
bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_id_string()==fid2.get_id_string());
}

bool operator< (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_id_string()<fid2.get_id_string());
}

} // namespace scream
