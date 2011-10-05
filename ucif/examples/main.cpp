#include <fstream>
#include <iostream>
#include <vector>
#include <ucif/parser.h>
#include <ucif/builder.h>

namespace ucif { namespace example {

struct my_array_wrapper : array_wrapper_base
{
  std::vector<std::string> array;

  my_array_wrapper()
  : array()
  {}

  virtual void push_back(std::string const& value)
  {
    array.push_back(value);
  }

  virtual std::string operator[](unsigned const& i)
  {
    return array[i];
  }

  virtual unsigned size()
  {
    return array.size();
  }
};

struct my_builder : builder_base
{
  virtual void start_save_frame(std::string const& save_frame_heading) {}
  virtual void end_save_frame() {}
  virtual void add_data_item(std::string const& tag, std::string const& value) {}
  virtual void add_loop(array_wrapper_base const& loop_headers,
                        array_wrapper_base const& values) {}
  virtual void add_data_block(std::string const& data_block_heading) {}
  virtual array_wrapper_base* new_array()
  {
    return new my_array_wrapper();
  }
};

}} // namespace ucif::example

int main (int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "Please provide a path to a CIF file." << std::endl;
    return 0;
  }
  std::string filename(argv[1]);
  std::string input_string;
  std::ifstream myfile(argv[1], std::ifstream::in);
  if (!myfile.is_open()) {
    std::cout << "Could not open file " << argv[1] << std::endl;
    return 0;
  }
  std::string tmp;
  while (getline(myfile, tmp)) {
    input_string += tmp;
    input_string += "\n";
  }
  myfile.close();

  ucif::example::my_builder builder;
  ucif::parser parsed(&builder, input_string, filename, /*strict=*/true);

  // Were there any lexing/parsing/tree walking errors?
  std::vector<std::string> lexer_errors =
    dynamic_cast<ucif::example::my_array_wrapper*>(parsed.lxr->errors)->array;
  std::vector<std::string> parser_errors =
    dynamic_cast<ucif::example::my_array_wrapper*>(parsed.psr->errors)->array;
  std::vector<std::string> tree_walker_errors;
  if (parsed.tree_psr != NULL) {
    tree_walker_errors =
      dynamic_cast<ucif::example::my_array_wrapper*>(parsed.tree_psr->errors)->array;
  }
  for (int i=0;i<lexer_errors.size();i++) {
    std::cout << lexer_errors[i] << std::endl;
  }
  for (int i=0;i<parser_errors.size();i++) {
    std::cout << parser_errors[i] << std::endl;
  }
  for (int i=0;i<tree_walker_errors.size();i++) {
    std::cout << tree_walker_errors[i] << std::endl;
  }
  if (lexer_errors.size() + parser_errors.size() + tree_walker_errors.size() == 0) {
    std::cout << "Congratulations! " << argv[1] <<
    " is a syntactically correct CIF file!" << std::endl;
  }

  return 0;
}
