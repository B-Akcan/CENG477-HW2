#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

double f_(double x, double y, double x_0, double y_0, double x_1, double y_1) {
    return (x * (y_0 - y_1)) + (y * (x_1 - x_0)) + (x_0 * y_1) - (y_0 * x_1);
}

double calculateYMin(Vec4 v0, Vec4 v1, Vec4 v2) {
	return std::min({v0.y, v1.y, v2.y});
}

double calculateXMin(Vec4 v0, Vec4 v1, Vec4 v2) {
	return std::min({v0.x, v1.x, v2.x});
}

double calculateYMax(Vec4 v0, Vec4 v1, Vec4 v2) {
	return std::max({v0.y, v1.y, v2.y});
}

double calculateXMax(Vec4 v0, Vec4 v1, Vec4 v2) {
	return std::max({v0.x, v1.x, v2.x});
}

int roundDoubleToInt(double val) {
	return (int)(val + 0.5);
}

Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	if(depthBuffer.empty()){
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<double> rowOfBuffer;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfBuffer.push_back(1.1);
			}

			this->depthBuffer.push_back(rowOfBuffer);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
				this->depthBuffer[i][j] = __DBL_MAX__;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	try {
		fout.open(camera->outputFilename);
	} catch (const std::string& ex) {
		std::cout << ex << std::endl;
	}
	

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function
	Matrix4 M_cam = calculateCameraTransformationMatrix(camera);
	Matrix4 M_proj = calculateProjectionTransformationMatrix(camera);
	Matrix4 M_vp = calculateViewportTransformationMatrix(camera);
	

	for (Mesh *mesh : this->meshes) {
		Matrix4 M_modeling = calculateModelingTransformationMatrix(camera, mesh);
		Matrix4 M_proj_cam_modeling = multiplyMatrixWithMatrix(M_proj, multiplyMatrixWithMatrix(M_cam, M_modeling));

		for (int j = 0; j < mesh->numberOfTriangles; j++) {
			Triangle triangle = mesh->triangles[j];
			Vec3 *v0 = this->vertices[triangle.vertexIds[0]-1];
			Vec3 *v1 = this->vertices[triangle.vertexIds[1]-1];
			Vec3 *v2 = this->vertices[triangle.vertexIds[2]-1];
			Color v0_color = *(this->colorsOfVertices[v0->colorId-1]);
			Color v1_color = *(this->colorsOfVertices[v1->colorId-1]);
			Color v2_color = *(this->colorsOfVertices[v2->colorId-1]);

			Vec4 v0_homo = {v0->x, v0->y, v0->z, 1, v0->colorId};
			Vec4 v1_homo = {v1->x, v1->y, v1->z, 1, v1->colorId};
			Vec4 v2_homo = {v2->x, v2->y, v2->z, 1, v2->colorId};

			Vec4 v0_transformed = multiplyMatrixWithVec4(M_proj_cam_modeling, v0_homo);
			Vec4 v1_transformed = multiplyMatrixWithVec4(M_proj_cam_modeling, v1_homo);
			Vec4 v2_transformed = multiplyMatrixWithVec4(M_proj_cam_modeling, v2_homo);

			if (this->cullingEnabled && backfaceCulling(v0_transformed, v1_transformed, v2_transformed) < 0) {
				continue;
			}

			// Perspective division (4th element becomes 1, others are /t)
			Vec4 v0_perspective_divided = {v0_transformed.x / v0_transformed.t, v0_transformed.y / v0_transformed.t, v0_transformed.z / v0_transformed.t, 1, v0_transformed.colorId};
			Vec4 v1_perspective_divided = {v1_transformed.x / v1_transformed.t, v1_transformed.y / v1_transformed.t, v1_transformed.z / v1_transformed.t, 1, v1_transformed.colorId};
			Vec4 v2_perspective_divided = {v2_transformed.x / v2_transformed.t, v2_transformed.y / v2_transformed.t, v2_transformed.z / v2_transformed.t, 1, v2_transformed.colorId};

			if (mesh->type == 0) {
				//prevent double Vec4 and Color modification from multiple LiangBarsky calling
				Vec4 v0_copy = v0_perspective_divided;
				Vec4 v1_copy = v1_perspective_divided;
				Color v0_color_copy = v0_color;
				Color v1_color_copy = v1_color;
				if(liangBarsky(v0_copy, v1_copy, v0_color_copy, v1_color_copy)){
					Vec4 v0_viewport = multiplyMatrixWithVec4(M_vp, v0_copy);
					Vec4 v1_viewport = multiplyMatrixWithVec4(M_vp, v1_copy);
					
					rasterizeLine(this->image, v0_viewport, v1_viewport, v0_color_copy, v1_color_copy);
				}
				Vec4 v1_copy_second = v1_perspective_divided;
				Vec4 v2_copy = v2_perspective_divided;
				Color v1_color_copy_second = v1_color;
				Color v2_color_copy = v2_color;
				if(liangBarsky(v1_copy_second, v2_copy, v1_color_copy_second, v2_color_copy)){
					Vec4 v1_viewport = multiplyMatrixWithVec4(M_vp, v1_copy_second);
					Vec4 v2_viewport = multiplyMatrixWithVec4(M_vp, v2_copy);
					
					rasterizeLine(this->image, v1_viewport, v2_viewport, v1_color_copy_second, v2_color_copy);
				}
				Vec4 v2_copy_second = v2_perspective_divided;
				Vec4 v0_copy_second = v0_perspective_divided;
				Color v2_color_copy_second = v2_color;
				Color v0_color_copy_second = v0_color;
				if(liangBarsky(v2_copy_second, v0_copy_second, v2_color_copy_second, v0_color_copy_second)){
					Vec4 v2_viewport = multiplyMatrixWithVec4(M_vp, v2_copy_second);
					Vec4 v0_viewport = multiplyMatrixWithVec4(M_vp, v0_copy_second);
					
					rasterizeLine(this->image, v2_viewport, v0_viewport, v2_color_copy_second, v0_color_copy_second);
				}
			}
			else {
				Vec4 v0_viewport = multiplyMatrixWithVec4(M_vp, v0_perspective_divided);
				Vec4 v1_viewport = multiplyMatrixWithVec4(M_vp, v1_perspective_divided);
				Vec4 v2_viewport = multiplyMatrixWithVec4(M_vp, v2_perspective_divided);

				rasterizeTriangle(this->image, v0_viewport, v1_viewport, v2_viewport, v0_color, v1_color, v2_color);
			}
		}
	}
}

Matrix4 Scene::calculateCameraTransformationMatrix(Camera *camera) {
    double T[4][4] = {{1, 0, 0, -(camera->position.x)},
                        {0, 1, 0, -(camera->position.y)},
                        {0, 0, 1, -(camera->position.z)},
                        {0, 0, 0, 1}};
    double R[4][4] = {{camera->u.x, camera->u.y, camera->u.z, 0},
                      {camera->v.x, camera->v.y, camera->v.z, 0},
                      {camera->w.x, camera->w.y, camera->w.z, 0},
                      {0, 0, 0, 1}};

	Matrix4 T_mat = Matrix4(T);
	Matrix4 R_mat = Matrix4(R);
    return multiplyMatrixWithMatrix(R_mat, T_mat);
}

Matrix4 Scene::calculateProjectionTransformationMatrix(Camera *camera) {
	if (camera->projectionType == 0) { // orthographic
		double M[4][4] = {{2.0 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
							{0, 2.0 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
							{0, 0, -2.0 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
							{0, 0, 0, 1}};
		return Matrix4(M);
	} else if (camera->projectionType == 1) { // perspective
		double M[4][4] = {{2.0 * camera->near / (camera->right - camera->left), 0, (camera->right + camera->left) / (camera->right - camera->left), 0},
							{0, (2.0 * camera->near) / (camera->top - camera->bottom), (camera->top + camera->bottom) / (camera->top - camera->bottom), 0},
							{0, 0, -(camera->far + camera->near) / (camera->far - camera->near), (-2.0 * camera->far * camera->near) / (camera->far - camera->near)},
							{0, 0, -1, 0}};
		return Matrix4(M);
	}

	return Matrix4();
}

Matrix4 Scene::calculateViewportTransformationMatrix(Camera *camera) {
    double M_viewport[4][4] = {{camera->horRes/2.0, 0, 0, (camera->horRes-1)/2.0},
                               {0, camera->verRes/2.0, 0, (camera->verRes-1)/2.0},
                               {0, 0, 0.5, 0.5},
                               {0, 0, 0, 1}};
    return Matrix4(M_viewport);
}

Matrix4 Scene::calculateModelingTransformationMatrix(Camera *camera, Mesh *mesh) {
	Matrix4 M_modeling = getIdentityMatrix();
	for (int i = 0; i < mesh->numberOfTransformations; i++) {
		char transformationType = mesh->transformationTypes[i];
		int transformationId = mesh->transformationIds[i] -1; // since transformations start from 1
		

		if (transformationType == 't') { // translation
			double tx = this->translations[transformationId]->tx;
			double ty = this->translations[transformationId]->ty;
			double tz = this->translations[transformationId]->tz;

			double M_translation[4][4] = {{1, 0, 0, tx},
									{0, 1, 0, ty},
									{0, 0, 1, tz},
									{0, 0, 0, 1}};

			M_modeling = multiplyMatrixWithMatrix(M_translation, M_modeling);
		} else if (transformationType == 's') { // scaling
			double sx = this->scalings[transformationId]->sx;
			double sy = this->scalings[transformationId]->sy;
			double sz = this->scalings[transformationId]->sz;

			double M_scaling[4][4] = {{sx, 0, 0, 0},
										{0, sy, 0, 0},
										{0, 0, sz, 0},
										{0, 0, 0, 1}};
									
			M_modeling = multiplyMatrixWithMatrix(M_scaling, M_modeling);
		} else if (transformationType == 'r') { // rotation
			double angle = this->rotations[transformationId]->angle;
			double ux = this->rotations[transformationId]->ux;
			double uy = this->rotations[transformationId]->uy;
			double uz = this->rotations[transformationId]->uz;
            
            Vec3 u = Vec3(ux, uy, uz, -1), v, w;
			double abs_ux = abs(ux);
			double abs_uy = abs(uy);
			double abs_uz = abs(uz);
            double minimumComponent = std::min(std::min(abs_ux, abs_uy), abs_uz);
            if (minimumComponent == abs_ux)
                v = Vec3(0, -uz, uy, -1);
            else if (minimumComponent == abs_uy)
                v = Vec3(-uz, 0, ux, -1);
            else if (minimumComponent == abs_uz)
                v = Vec3(-uy, ux, 0, -1);
            w = crossProductVec3(u, v);
            
            v = normalizeVec3(v);
            w = normalizeVec3(w);

            double M[4][4] = {{u.x, u.y, u.z, 0},
								{v.x, v.y, v.z, 0},
								{w.x, w.y, w.z, 0},
								{0, 0, 0, 1}};
            double M_inverse[4][4] = {{u.x, v.x, w.x, 0},
										{u.y, v.y, w.y, 0},
										{u.z, v.z, w.z, 0},
										{0, 0, 0, 1}};
            
			double radian_angle = angle * M_PI/180;
            double R_x_theta[4][4] = {{1, 0, 0, 0},
                                    {0, cos(radian_angle), -sin(radian_angle), 0},
                                    {0, sin(radian_angle), cos(radian_angle), 0},
                                    {0, 0, 0, 1}};
            Matrix4 temp = multiplyMatrixWithMatrix(R_x_theta, M);
            temp = multiplyMatrixWithMatrix(M_inverse, temp);
            M_modeling = multiplyMatrixWithMatrix(temp, M_modeling);
		}
	}

	return M_modeling;
}

double Scene::backfaceCulling(Vec4 v0_transformed, Vec4 v1_transformed, Vec4 v2_transformed) {
	Vec3 v1_minus_v0 = subtractVec3({v1_transformed.x, v1_transformed.y, v1_transformed.z}, {v0_transformed.x, v0_transformed.y, v0_transformed.z});
	Vec3 v2_minus_v0 = subtractVec3({v2_transformed.x, v2_transformed.y, v2_transformed.z}, {v0_transformed.x, v0_transformed.y, v0_transformed.z});
	Vec3 normal = crossProductVec3(v1_minus_v0, v2_minus_v0);
	Vec3 v = {v0_transformed.x, v0_transformed.y, v0_transformed.z}; // all vertices are act as camera eye moved to (0,0,0)
	double n_dot_v = dotProductVec3(normal, v);

	return n_dot_v;
}

bool Scene::visible(double den, double num, double &t_e, double &t_l) {
	if (den > 0) { // potentially entering
		double t = num / den;
		if (t > t_l)
			return false;
		if (t > t_e)
			t_e = t;
	}
	else if (den < 0) { // potentially leaving
		double t = num / den;
		if (t < t_e)
			return false;
		if (t < t_l)
			t_l = t;
	}
	else if (num > 0) // line parallel to edge
		return false;

	return true;
}

//in place modification of vertices and colors
bool Scene::liangBarsky(Vec4 &v0, Vec4 &v1, Color& v0_color, Color& v1_color) {
	double t_e = 0, t_l = 1;
	bool isVisible = false;
    
	double dx = v1.x - v0.x;
	double dy = v1.y - v0.y;
	double dz = v1.z - v0.z;
    
    Color dc = v1_color - v0_color; //color change rate along the line
    double x_min = -1, y_min = -1, z_min = -1; //min normalized box borders

    double x_max = 1, y_max = 1, z_max = 1;
    if (visible(dx, x_min-v0.x, t_e, t_l) && visible(-dx, v0.x-x_max, t_e, t_l)
        && visible(dy, y_min-v0.y, t_e, t_l) && visible(-dy, v0.y-y_max, t_e, t_l)
        && visible(dz, z_min-v0.z, t_e, t_l) && visible(-dz, v0.z-z_max, t_e, t_l)) {
        isVisible = true;
        /* Check if at least some part of the line is clipped */
        if (t_l < 1) {
            v1.x = v0.x + (dx * t_l);
            v1.y = v0.y + (dy * t_l);
            v1.z = v0.z + (dz * t_l);
            v1_color = v0_color + (dc * t_l);
        }
        if (t_e > 0) {
            v0.x = v0.x + (dx * t_e);
            v0.y = v0.y + (dy * t_e);
            v0.z = v0.z + (dz * t_e);
            v0_color = v0_color + (dc * t_l);
        }
    }


    return isVisible;
}

void Scene::rasterizeLine(vector<vector<Color>> &image, Vec4 v0, Vec4 v1, Color c0, Color c1) {
	double abs_dx = abs(v1.x - v0.x);
	double abs_dy = abs(v1.y - v0.y);

    if (abs_dy <= abs_dx) {
		int sign;

        if (v1.x < v0.x) {
            swap(v0, v1);
            swap(c0, c1);
        }
        if (v1.y < v0.y) {
            sign = -1;
        }
		else {
			sign = 1;
		}

        int y = v0.y;
        Color c = c0;
		Color dc = (c1 - c0) / (v1.x - v0.x);
        double d = (v0.y - v1.y) + (sign * (v1.x - v0.x) / 2);
		double dz = (v1.z - v0.z) / (v1.x - v0.x);
		double current_z = v0.z;

        for (int x = v0.x; x <= v1.x + EPSILON; x++) {
			if (depthBuffer[x][y] > current_z) {
            	image[x][y] = c.round();
				depthBuffer[x][y] = current_z;
			}

            if (sign * d < 0) {
                y += sign;
                d += (v0.y - v1.y) + (sign * (v1.x - v0.x));
            }
            else
                d += (v0.y - v1.y);

            c = c + dc;
			current_z = current_z + dz;
        }
    }
    else {
		int sign;

        if (v1.y < v0.y) {
            swap(v0, v1);
            swap(c0, c1);
        }
        if (v1.x < v0.x) {
            sign = -1;
        }
		else {
			sign = 1;
		}

        int x = v0.x;
        Color c = c0;
		Color dc = (c1 - c0) / (v1.y - v0.y);
        double d = (v1.x - v0.x) + (sign * (v0.y - v1.y) / 2);
		double dz = (v1.z - v0.z) / (v1.y - v0.y);
		double current_z = v0.z;

        for (int y = v0.y; y <= v1.y + EPSILON; y++) {
            if (depthBuffer[x][y] > current_z) {
            	image[x][y] = c.round();
				depthBuffer[x][y] = current_z;
			}

            if (sign * d > 0) {
                x += sign;
                d += (v1.x - v0.x) + (sign * (v0.y - v1.y));
            }
            else
                d += (v1.x - v0.x);
				
            c = c + dc;
			current_z = current_z + dz;
        }
    }
}

void Scene::rasterizeTriangle(vector<vector<Color>> &image, Vec4 v0, Vec4 v1, Vec4 v2, Color c0, Color c1, Color c2) {
	int ymin = (calculateYMin(v0, v1, v2));
	int xmin = (calculateXMin(v0, v1, v2));
	int ymax = (calculateYMax(v0, v1, v2));
	int xmax = (calculateXMax(v0, v1, v2));
	
	for (int y = ymin; y <= ymax; y++) {
		for (int x = xmin; x <= xmax; x++) {
			if(x < 0 || y < 0)
				continue;
			double alpha = f_(x, y, v1.x, v1.y, v2.x, v2.y) / f_(v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
			double beta = f_(x, y, v2.x, v2.y, v0.x, v0.y) / f_(v1.x, v1.y, v2.x, v2.y, v0.x, v0.y);
			double gama = f_(x, y, v0.x, v0.y, v1.x, v1.y) / f_(v2.x, v2.y, v0.x, v0.y, v1.x, v1.y);
			if (alpha >= 0 && beta >= 0 && gama >= 0 ) {
				Color c = c0 * alpha + c1 * beta + c2 * gama;
				double z_value = alpha * v0.z + beta * v1.z + gama * v2.z;
				double value = depthBuffer[x][y];
				if(value > z_value - EPSILON){
					image[x][y] = c.round();
					depthBuffer[x][y] = z_value;
				}
			}
		}
	}
}

