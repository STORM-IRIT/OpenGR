#include <Eigen/Core>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/time.h>
#include <pcl/console/print.h>
#include <pcl/features/fpfh_omp.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <string>
#include <functional> // std::function
#include <map>

#include <pcl/registration/super4pcs.h>

#include <gr/shared.h>
#include "../demo-utils.h"

// Types
typedef pcl::PointNormal PointNT;
typedef pcl::PointCloud<PointNT> PointCloudT;
typedef pcl::visualization::PointCloudColorHandlerCustom<PointNT> ColorHandlerT;

using namespace gr;

using loadfunc = std::function<int(const std::string &, pcl::PointCloud<PointNT> &)>;
const std::map<std::string, loadfunc> loaders =
{
    { "obj", pcl::io::loadOBJFile<PointNT> },
    { "ply", pcl::io::loadPLYFile<PointNT> },
    { "pcd", pcl::io::loadPCDFile<PointNT> }
};


// src: https://en.cppreference.com/w/cpp/string/byte/tolower
std::string str_tolower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); } // correct
                  );
    return s;
}

bool
load(const std::string& filename, PointCloudT::Ptr& pcloud){
    auto getFileExt = [] (const std::string& s) -> std::string {

       size_t i = s.rfind('.', s.length());
       return i != std::string::npos
                 ? str_tolower( s.substr(i+1, s.length() - i) )
                 : "";
    };

    std::string ext = getFileExt(filename);
    auto l = loaders.find(ext);
    if( l != loaders.end() ) {
        return (*l).second( filename, *pcloud ) >= 0;
    }
    pcl::console::print_error ("Unsupported file extension: %s\n", ext.c_str());
    return false;
}


// Align a rigid object to a scene with clutter and occlusions
int
main (int argc, char **argv)
{
  // Point clouds
  PointCloudT::Ptr object (new PointCloudT);
  PointCloudT::Ptr object_aligned (new PointCloudT);
  PointCloudT::Ptr scene (new PointCloudT);

  // Get input object and scene
  if (argc < 3)
  {
    pcl::console::print_error ("Syntax is: %s scene.obj object.obj [PARAMS]\n", argv[0]);
    Demo::printParameterList();
    return (-1);
  }

  std::string objPath {argv[2]};
  std::string scenePath {argv[1]};

  // Load object and scene
  pcl::console::print_highlight ("Loading point clouds...\n");
  if ( ! ( load(objPath, object) && load(scenePath, scene) ) )
  {
    pcl::console::print_error ("Error loading object/scene file!\n");
    return (-1);
  }

  if(int c = Demo::getArgs(argc, argv) != 0)
    {
      Demo::printUsage(argc, argv);
      exit(std::max(c,0));
    }

  pcl::Super4PCS<PointNT,PointNT> align;
  auto &options = align.getOptions();
  Demo::setOptionsFromArgs(options);

  // Perform alignment
  pcl::console::print_highlight ("Starting alignment...\n");
  align.setInputSource (object);
  align.setInputTarget (scene);

  {
    pcl::ScopeTime t("Alignment");
    align.align (*object_aligned);
  }

  if (align.hasConverged ())
  {
    // Print results
    printf ("\n");
    Eigen::Matrix4f transformation = align.getFinalTransformation ();
    pcl::console::print_info ("    | %6.3f %6.3f %6.3f | \n", transformation (0,0), transformation (0,1), transformation (0,2));
    pcl::console::print_info ("R = | %6.3f %6.3f %6.3f | \n", transformation (1,0), transformation (1,1), transformation (1,2));
    pcl::console::print_info ("    | %6.3f %6.3f %6.3f | \n", transformation (2,0), transformation (2,1), transformation (2,2));
    pcl::console::print_info ("\n");
    pcl::console::print_info ("t = < %0.3f, %0.3f, %0.3f >\n", transformation (0,3), transformation (1,3), transformation (2,3));
    pcl::console::print_info ("\n");

    // Show alignment
    pcl::visualization::PCLVisualizer visu("Alignment - Super4PCS");
    visu.addPointCloud (scene, ColorHandlerT (scene, 0.0, 255.0, 0.0), "scene");
    visu.addPointCloud (object_aligned, ColorHandlerT (object_aligned, 0.0, 0.0, 255.0), "object_aligned");

    pcl::console::print_highlight ("Saving registered cloud to %s ...\n", Demo::defaultPlyOutput.c_str());
    pcl::io::savePLYFileBinary<PointNT>(Demo::defaultPlyOutput, *object_aligned);

    visu.spin ();
  }
  else
  {
    pcl::console::print_error ("Alignment failed!\n");
    return (-1);
  }

  return (0);
}
