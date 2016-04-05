#include <iostream>
#include <getopt.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#ifdef USE_THREADS
#include <queue>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#endif

#include "./qm_interpolation.h"

using std::string;

static void usage()
{
  using std::cerr;
  using std::endl;
  cerr << endl;
  cerr << "About: energy force calculation, current only support waterbox" << endl;
  cerr << "Usage: qm [option]" << endl;
  cerr << "Options:" << endl;
  cerr << "  -d, --database  the database file(gz file)" << endl;
  cerr << "  -i, --input     input dir which contains .inp files" << endl;
  cerr << "  -o, --ouput     output file" << endl;
  cerr << "  -m, --model     current only support \"wtr\"" << endl;
#ifdef USE_THREADS
  cerr << "  -t, --threads   threads count" << endl;
#endif

  cerr << "Example:" << endl;
  cerr << "  qm -d database -m wtr -i random -o output" << endl;
#ifdef USE_THREADS
  cerr << "  qm -t 2 -d database -m wtr -i random -o output" << endl;
#endif

  exit(1);
}

#ifdef USE_THREADS

boost::mutex io_mutex;
boost::mutex queue_mutex;
std::queue<string> filenames;

void process(FILE*& fp_out, QM::QMInterpolation& interpolation)
{
  string filename;
  while (true)
  {
    {
      boost::mutex::scoped_lock lock(queue_mutex);
      if (filenames.empty())
      {
        break;
      }
      filename = filenames.front();
      filenames.pop();
    }
    const char* output = interpolation.process(filename).c_str();
    {
      boost::mutex::scoped_lock lock(io_mutex);
      fprintf(fp_out, "%s\n", output);
    }
  }
}

#endif

int main(int argc, char* argv[])
{

  string database_file;
  string input_dir;
  string output_file;
  string model;

#ifdef USE_THREADS
  int t = 1;
#endif

  int c;
  static struct option loptions[] =
  {
    {"help", no_argument, NULL, 'h'},
    {"database", required_argument, NULL, 'd'},
    {"input", required_argument, NULL, 'i'},
    {"ouput", required_argument, NULL, 'o'},
    {"model", required_argument, NULL, 'm'},
#ifdef USE_THREADS
    {"threads", optional_argument, NULL, 't'},
#endif
    {NULL, 0, NULL, 0}
  };

#ifdef USE_THREADS
  while ((c = getopt_long(argc, argv, "h?d:i:o:m:t:?", loptions, NULL)) != -1)
#else
  while ((c = getopt_long(argc, argv, "h?d:i:o:m:", loptions, NULL)) != -1)
#endif
  {
    switch(c)
    {
    case 'd':
      database_file = optarg;
      break;
    case 'i':
      input_dir = optarg;
      break;
    case 'o':
      output_file = optarg;
      break;
    case 'm':
      model = optarg;
      break;
#ifdef USE_THREADS
    case 't':
      t = boost::lexical_cast<int>(optarg);
      break;
#endif
    default:
      usage();
    }
  }
  if (optind < 9) // 4 args and qm
  {
    std::cerr << "Error: you input less argv, see:" << std::endl;
    usage();
  }

  if (model != "wtr")
  {
    std::cerr << "current only support waterbox, you input should be wtr" << std::endl;
    exit(1);
  }
  database::EnergeForceDatabase energe(database_file.c_str(), model.c_str());
  energe.read_file();
  QM::QMInterpolation interpolation(model, energe);

  boost::filesystem::path path(input_dir);
  boost::filesystem::directory_iterator end_iter;

  FILE* fp_out = fopen(output_file.c_str(), "w");
  for (boost::filesystem::directory_iterator iter(path); iter!=end_iter; ++iter)
  {
    if (boost::filesystem::is_regular_file(iter->status()) && boost::ends_with(iter->path().filename().string(), ".inp"))
    {
#ifdef USE_THREADS
       filenames.push(iter->path().string());
#else
       fprintf(fp_out, "%s\n", interpolation.process(iter->path().string()).c_str());
#endif
    }
  }

#ifdef USE_THREADS
  if (t <= 0)
  {
    std::cerr << "Error: theads should be a postive value" << std::endl;
    exit(1);
  }
  boost::thread_group threads;
  for(int i=0; i<t; i++)
  {
    threads.create_thread(boost::bind(process, fp_out, interpolation));
  }

  threads.join_all();
#endif

  fclose(fp_out);
  return 0;
}
