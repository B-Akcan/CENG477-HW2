#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double>>depthBuffer;//between -1 and 1 doubles
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);
	void forwardRenderingPipeline(Camera *camera);
	Matrix4 calculateCameraTransformationMatrix(Camera * camera);
	Matrix4 calculateProjectionTransformationMatrix(Camera *camera);
	Matrix4 calculateViewportTransformationMatrix(Camera * camera);
	Matrix4 calculateModelingTransformationMatrix(Camera *camera, Mesh *mesh);
	double backfaceCulling(Vec4 v0_transformed, Vec4 v1_transformed, Vec4 v2_transformed);
	bool visible(double den, double num, double &t_e, double &t_l);
	bool liangBarsky(Vec4 &v0, Vec4 &v1, Color& v0_color, Color& v1_color);
	void rasterizeLine(std::vector<std::vector<Color>> &image, Vec4 v0, Vec4 v1, Color c0, Color c1);
	void rasterizeTriangle(std::vector<std::vector<Color>> &image, Vec4 v0, Vec4 v1, Vec4 v2, Color c0, Color c1, Color c2);
};

#endif
