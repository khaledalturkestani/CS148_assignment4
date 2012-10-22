/****************************************************************************
 * THE GRAND METAMORPHOSIS
 * CS148 Assignment #4 - Fall 2010, Stanford University
 ****************************************************************************/

#include "st.h"
#include "stglut.h"
#include "parseConfig.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// --------------------------------------------------------------------------
// Structure to contain an image feature for a morph. A feature is a directed
// line segment from P to Q, with coordinates in pixel units relative to the
// lower-left corner of the image.
// --------------------------------------------------------------------------

struct Feature
{
    STPoint2 P, Q;
    Feature(const STPoint2 &p, const STPoint2 &q) : P(p), Q(q) { }
};

// --------------------------------------------------------------------------
// Constants, a few global variables, and function prototypes
// --------------------------------------------------------------------------

const int kWindowWidth  = 512;
const int kWindowHeight = 512;
const int kFrames       = 30;   // number of frames to generate

STImage *gDisplayedImage = 0;   // an image to display (for testing/debugging)

std::vector<Feature> gSourceFeatures;   // feature set on source image
std::vector<Feature> gTargetFeatures;   // corresponding features on target

// Copies an image into the global image for display
void DisplayImage(STImage *image);

// --------------------------------------------------------------------------
// CS148 TODO: Implement the functions below to compute the morph
// --------------------------------------------------------------------------

/**
 * Compute a linear blend of the pixel colors in two provided images according
 * to a parameter t.
 */
STImage *BlendImages(STImage *image1, STImage *image2, float t)
{
    STColor4ub image1Color;
    STColor4ub image2Color;
    STImage *result = new STImage(*image1);
    for (int i = 0; i < kWindowWidth; i++) {
        for (int j = 0; j < kWindowHeight; j++) {
            image1Color = image1->GetPixel(i, j);
            image2Color = image2->GetPixel(i, j);
            int newR = image1Color.r + t*(image2Color.r-image1Color.r);
            int newG = image1Color.g + t*(image2Color.g-image1Color.g);
            int newB = image1Color.b + t*(image2Color.b-image1Color.b);
            result->SetPixel(i, j, STColor4ub(newR, newG, newB));
        }
    }    
    return result;
}

STColor4ub interpolateColor(STColor4ub pixel1, STColor4ub pixel2, float t)
{
    int newR = pixel1.r + t*(pixel2.r-pixel1.r);
    int newG = pixel1.g + t*(pixel2.g-pixel1.g);
    int newB = pixel1.b + t*(pixel2.b-pixel1.b);
    return STColor4ub(newR, newG, newB);
}

STColor4ub bilinear(STImage *image, STVector2 point)
{
    float s = point.x - (int)point.x;
    float t = point.y - (int)point.y;
    STColor4ub v0 = image->GetPixel((int)point.x, (int)point.y);
    STColor4ub v1 = image->GetPixel((int)point.x+1, (int)point.y);
    STColor4ub v2 = image->GetPixel((int)point.x, (int)point.y+1);
    STColor4ub v3 = image->GetPixel((int)point.x+1, (int)point.y+1);
    STColor4ub v01 = interpolateColor(v0, v1, s);
    STColor4ub v23 = interpolateColor(v2, v3, s);
    STColor4ub v = interpolateColor(v01, v23, t);
    return v;
}

STVector2 Perpendicular(STVector2 vec)
{
    STVector2 perpendicular = STVector2(vec.y, (-1)*vec.x);
    return perpendicular;
}

float shortestDistance(STVector2 P, STVector2 Q, STVector2 X, float u, float v)
{
    float distance = 0;
    if (u >= 0.0 && u <= 1.0) {
        distance = abs(v);
    } else if (u < 0.0) {
        distance = (P-X).Length();
    } else {
        distance = (Q-X).Length();
    }
    return distance;
}
/**
 * Compute a field morph on an image using two sets of corresponding features
 * according to a parameter t.  Arguments a, b, and p are weighting parameters
 * for the field morph, as described in Beier & Nelly 1992, section 3.
 */
STImage *FieldMorph(STImage *image,
                    const std::vector<Feature> &sourceFeatures,
                    const std::vector<Feature> &targetFeatures,
                    float t, float a, float b, float p)
{
    STImage *result = new STImage(kWindowWidth, kWindowHeight);
    for (int i = 0; i < kWindowWidth; i++) {
        for (int j = 0; j < kWindowHeight; j++) {
            STVector2 DSUM = STVector2(0, 0);
            float weightSum = 0;
            STPoint2 Xpoint = STPoint2(i, j);
            STVector2 X = STVector2(Xpoint);;
            for (int k = 0; k < sourceFeatures.size(); k++) {
                STVector2 P = STVector2(targetFeatures[k].P);
                STVector2 Q = STVector2(targetFeatures[k].Q);
                STVector2 Pprime = STVector2(sourceFeatures[k].P);
                STVector2 Qprime = STVector2(sourceFeatures[k].Q);
                float u = STVector2::Dot(X-P, Q-P)/(Q-P).LengthSq();
                float v = STVector2::Dot(X-P, Perpendicular(Q-P))/(Q-P).Length();
                STVector2 Xp = Pprime + u*(Qprime-Pprime) + (v*Perpendicular(Qprime-Pprime))/(Qprime-Pprime).Length();
                STVector2 D = Xp - X;
                float dist = shortestDistance(P, Q, X, u, v);
                float weight = powf(powf((Q-P).Length(), p)/(a+dist), b);
                DSUM += D*weight;
                weightSum += weight;
            }
            STVector2 Xprime = X + DSUM/weightSum;

            // Adjust for out of bound coordinates
            if (Xprime.x >= kWindowWidth) Xprime.x = kWindowWidth-1;
            if (Xprime.x < 0) Xprime.x = 0;
            if (Xprime.y >= kWindowHeight) Xprime.y = kWindowHeight-1;
            if (Xprime.y < 0) Xprime.y = 0;
            STVector2 interpolatedVec = STVector2::Lerp(X, Xprime, t);

            if (Xprime.x >= 0 && Xprime.x < kWindowWidth-1 && Xprime.y >= 0 && Xprime.y < kWindowHeight-1)
            {
//                result->SetPixel(X.x, X.y, image->GetPixel(interpolatedVec.x, interpolatedVec.y));
                result->SetPixel(X.x, X.y, bilinear(image, interpolatedVec));
            } else {
                result->SetPixel(X.x, X.y, image->GetPixel(interpolatedVec.x, interpolatedVec.y));
            }
        }
    }
    return result;
}

/**
 * Compute a morph between two images by first distorting each toward the
 * other, then combining the results with a blend operation.
 */
STImage *MorphImages(STImage *sourceImage, const std::vector<Feature> &sourceFeatures,
                     STImage *targetImage, const std::vector<Feature> &targetFeatures,
                     float t, float a, float b, float p)
{
    STImage *src = FieldMorph(sourceImage, sourceFeatures, targetFeatures, t, a, b, p);
    STImage *dest = FieldMorph(targetImage, sourceFeatures, targetFeatures, 1.0-t, a, b, p);
    return BlendImages(src, dest, t);
}

/**
 * Compute a morph through time by generating appropriate values of t and
 * repeatedly calling MorphImages(). Saves the image sequence to disk.
 */
void GenerateMorphFrames(STImage *sourceImage, const std::vector<Feature> &sourceFeatures,
                         STImage *targetImage, const std::vector<Feature> &targetFeatures,
                         float a, float b, float p)
{
    // iterate and generate each required frame
    float t[] = {0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.18, 0.21, 0.26, 0.30, 0.36, 0.40, 0.46, 0.50,
                 0.55, 0.60, 0.65, 0.70, 0.73, 0.79, 0.81, 0.87, 0.90, 0.92, 0.94, 0.98, 0.99, 1.00};
    for (int i = 0; i < kFrames; ++i)
    {
        std::cout << "Metamorphosizing frame #" << i << "...";
        
        STImage *result = new STImage(*MorphImages(sourceImage, sourceFeatures, targetImage, targetFeatures, t[i], a, b, p));
        
        // generate a file name to save
        std::ostringstream oss;
        oss << "frame" << std::setw(3) << std::setfill('0') << i << ".png";

        // write and deallocate the morphed image
        if (result) {
            result->Save(oss.str());
            delete result;
        }

        std::cout << " done." << std::endl;
    }
}

// --------------------------------------------------------------------------
// Utility and support code below that you do not need to modify
// --------------------------------------------------------------------------

/**
 * Copies an image into the global image for display
 */
void DisplayImage(STImage *image)
{
    // clean up the previous image
    if (gDisplayedImage) {
        delete gDisplayedImage;
        gDisplayedImage = 0;
    }

    // allocate a new image and copy it over
    if (image) {
        gDisplayedImage = new STImage(image->GetWidth(), image->GetHeight());
        size_t bytes = image->GetWidth() * image->GetHeight() * sizeof(STImage::Pixel);
        memcpy(gDisplayedImage->GetPixels(), image->GetPixels(), bytes);
    }
}

/**
 * Display callback function draws a single image to help debug
 */
void DisplayCallback()
{
    glClearColor(.2f, 2.f, 2.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (gDisplayedImage)
        gDisplayedImage->Draw();

    glutSwapBuffers();
}

/**
 * Window resize callback function
 */
void ReshapeCallback(int w, int h)
{
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h);
}

/**
 * Keyboard callback function
 */
void KeyboardCallback(unsigned char key, int x, int y)
{
    switch (key)
    {
        // exit program on escape press
        case 27:
            exit(0);
            break;
        // save the currently displayed image if S is pressed
        case 's':
        case 'S':
            if (gDisplayedImage)
                gDisplayedImage->Save("screenshot.png");
            break;
        default:
            break;
    }
}

/**
 * This function is called by the parsing functions to populate the feature sets
 */
void AddFeatureCallback(STPoint2 p, STPoint2 q, ImageChoice image)
{
    if (image == IMAGE_1 || image == BOTH_IMAGES)
        gSourceFeatures.push_back(Feature(p, q));
    if (image == IMAGE_2 || image == BOTH_IMAGES)
        gTargetFeatures.push_back(Feature(p, q));
}

/**
 * Program entry point
 */
int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
    glutInitWindowPosition(20, 20);
    glutInitWindowSize(kWindowWidth, kWindowHeight);
    glutCreateWindow("Metamorphosis: CS148 Assignment 4");

    glutDisplayFunc(DisplayCallback);
    glutReshapeFunc(ReshapeCallback);
    glutKeyboardFunc(KeyboardCallback);

    //
    // load the configuration from config.txt, or other file as specified
    //
    std::string configFile = "config.txt";
    if (argc > 1) configFile = argv[1];

    char sourceName[64], targetName[64];
    char saveName[64], loadName[64];
    STImage *sourceImage, *targetImage;
    parseConfigFile(configFile.c_str(),
                    sourceName, targetName,
                    saveName, loadName,
                    &sourceImage, &targetImage);
    delete sourceImage;
    delete targetImage;

    //
    // load the features from the saved features file
    //
    loadLineEditorFile(loadName, AddFeatureCallback,
                       sourceName, targetName,
                       &sourceImage, &targetImage);

    //
    // run the full morphing algorithm before going into the main loop to
    // display an image
    //

    // these weighting parameters (Beier & Nelly 1992) can be changed if desired
    const float a = 0.5f, b = 1.0f, p = 0.2f;
    
    GenerateMorphFrames(sourceImage, gSourceFeatures,
                        targetImage, gTargetFeatures,
                        a, b, p);
     

    //
    // display a test or debug image here if desired
    // (note: comment this out if you call DisplayImage from elsewhere)
    //
    //STImage *result = sourceImage;

    // use this to test your image blending:
    //STImage *result = BlendImages(sourceImage, targetImage, 0.5f);
    

    // use this to test your field morph
    
    /*
    STImage *result = FieldMorph(sourceImage, gSourceFeatures, gTargetFeatures,
                                 0.0f, a, b, p);
    */
    // use this to test your image morphing
    /*
    STImage *result = MorphImages(sourceImage, gSourceFeatures,
                                  targetImage, gTargetFeatures,
                                  0.5f, a, b, p);
    */
    
    //DisplayImage(result);

    // enter the GLUT main loop
    glutMainLoop();

    return 0;
}

// --------------------------------------------------------------------------
