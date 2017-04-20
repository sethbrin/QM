#include <iostream>
#include <getopt.h>

#ifdef USE_THREADS
#include <queue>
#include <thread>
#include <mutex>
#endif

using std::string;

static void usage()
{
  using std::cerr;
  using std::endl;
  cerr << endl;
  cerr << "About: energy force calculation, current only support waterbox"
      << endl;
  cerr << "Usage: qm [option]" << endl;
  cerr << "Options:" << endl;
  cerr << "  -d, --database  the database file(gz file)" << endl;
  cerr << "  -i, --input     input dir which contains .inp files" << endl;
  cerr << "  -e, --ouput_energy     output file" << endl;
  cerr << "  -f, --ouput_force     output file" << endl;
  cerr << "  -t, --ouput_torque     output file" << endl;
#ifdef USE_THREADS
  cerr << "  -t, --threads   threads count" << endl;
#endif

  cerr << "Example:" << endl;
  cerr << "  qm -d database -i random -o output" << endl;
#ifdef USE_THREADS
  cerr << "  qm -t 2 -d database -i random -o output" << endl;
#endif

  exit(1);
}

#ifdef USE_THREADS

std::mutex io_mutex;
std::mutex queue_mutex;
std::queue<string> filenames;

//void process(FILE*& fp_out, QM::QMInterpolation& interpolation)
//{
//  string filename;
//  while (true)
//  {
//    {
//      boost::mutex::scoped_lock lock(queue_mutex);
//      if (filenames.empty())
//      {
//        break;
//      }
//      filename = filenames.front();
//      filenames.pop();
//    }
//    const char* output = interpolation.process(filename).c_str();
//    {
//      boost::mutex::scoped_lock lock(io_mutex);
//      fprintf(fp_out, "%s\n", output);
//    }
//  }
//}

#endif

int main(int argc, char* argv[])
{

  string database_file;
  string input_dir;
  string output_energy_file;
  string output_force_file;
  string output_torque_file;
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
    {"ouput_energy", required_argument, NULL, 'e'},
    {"ouput_force", required_argument, NULL, 'f'},
    {"ouput_torque", required_argument, NULL, 't'},
    {"ouput", required_argument, NULL, 'o'},
#ifdef USE_THREADS
    {"threads", optional_argument, NULL, '@'},
#endif
    {NULL, 0, NULL, 0}
  };

#ifdef USE_THREADS
  while ((c = getopt_long(argc, argv, "h?d:i:e:f:@:t:?", loptions, NULL)) != -1)
#else
  while ((c = getopt_long(argc, argv, "h?d:i:e:f:t:", loptions, NULL)) != -1)
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
    case 'e':
      output_energy_file = optarg;
      break;
    case 'f':
      output_force_file = optarg;
      break;
    case 't':
      output_torque_file = optarg;
      break;
    case 'm':
      model = optarg;
      break;
#ifdef USE_THREADS
    case '@':
      t = atoi(optarg);
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
    std::cerr << "current only support waterbox, you input should be wtr"
        << std::endl;
    exit(1);
  }
  //database::EnergeForceDatabase energe(database_file.c_str(), model.c_str());
  //energe.read_file();
  //QM::QMInterpolation interpolation(model, energe);

//  boost::filesystem::path path(input_dir);
//  boost::filesystem::directory_iterator end_iter;
//
//  FILE* fp_energy_out = fopen(output_energy_file.c_str(), "w");
//  FILE* fp_force_out = fopen(output_force_file.c_str(), "w");
//  FILE* fp_torque_out = fopen(output_torque_file.c_str(), "w");
//  for (boost::filesystem::directory_iterator iter(path);
//          iter!=end_iter; ++iter)
//  {
//    if (boost::filesystem::is_regular_file(iter->status()) &&
//            boost::ends_with(iter->path().filename().string(), ".inp"))
//    {
//#ifdef USE_THREADS
//       filenames.push(iter->path().string());
//#else
//       auto res = interpolation.process(iter->path().string());
//       fprintf(fp_energy_out, "%s\n",
//               res[0].c_str());
//       fprintf(fp_force_out, "%s\n",
//               res[1].c_str());
//       fprintf(fp_torque_out, "%s\n",
//               res[2].c_str());
//#endif
//    }
//  }
//
//#ifdef USE_THREADS
//  if (t <= 0)
//  {
//    std::cerr << "Error: theads should be a postive value" << std::endl;
//    exit(1);
//  }
//  boost::thread_group threads;
//  for(int i=0; i<t; i++)
//  {
//    threads.create_thread(boost::bind(process, fp_out, interpolation));
//  }
//
//  threads.join_all();
//#endif
//
//  fclose(fp_energy_out);
//  fclose(fp_force_out);
//  fclose(fp_torque_out);
  return 0;
}
