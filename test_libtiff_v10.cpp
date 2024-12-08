#include <tiffio.h>
#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#define u32 std::uint32_t
#define u16 std::uint16_t
#define u8 std::uint8_t
using namespace std;

/* 
********************************
Differences between v09 and v10:
********************************
thresholding

*/

// Function to read the TIFF image into a buffer
void readTiffImage(const char* filename, vector<u8>& buffer, u32& width, u32& height, u16& samples_per_pixel) {
    TIFF* tiff = TIFFOpen(filename, "r");
    if (tiff == nullptr) {
        cerr << "Could not open TIFF file!" << endl;
        return;
    }

    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

    buffer.resize(width * height * samples_per_pixel);
    for (u32 row = 0; row < height; row++) {
        TIFFReadScanline(tiff, &buffer[row * width * samples_per_pixel], row);
    }

    TIFFClose(tiff);
}


// Daubechies-4 low-pass and high-pass filter coefficients
const float h0 = (1 + sqrt(3)) / (4 * sqrt(2));
const float h1 = (3 + sqrt(3)) / (4 * sqrt(2));
const float h2 = (3 - sqrt(3)) / (4 * sqrt(2));
const float h3 = (1 - sqrt(3)) / (4 * sqrt(2));
const float g0 = -h3;
const float g1 = h2;
const float g2 = -h1;
const float g3 = h0;

// Coefficients for the inverse Daubechies-4 transform:
const float ih0 = h0;
const float ih1 = -h1;
const float ih2 = h2;
const float ih3 = -h3;

const float ig0 = g0;
const float ig1 = -g1;
const float ig2 = g2;
const float ig3 = -g3;

// Function to perform 1D Daubechies-4 transform on an array
void daubechies1D(vector<float>& data) {
    int N = data.size();
    vector<float> temp(N);

    for (int i = 0; i < N / 2; i++) {
        temp[i] = h0 * data[2*i] + h1 * data[2*i+1] + h2 * data[2*i+2] + h3 * data[2*i+3];
        temp[i + N/2] = g0 * data[2*i] + g1 * data[2*i+1] + g2 * data[2*i+2] + g3 * data[2*i+3];
    }
    for (int i = 0; i < N; i++) {
        data[i] = temp[i];
    }
}

// 2D Daubechies-4 Transform: Apply 1D transform on rows and columns
void daubechies2D(vector<vector<float>>& image) {
    int rows = image.size();
    int cols = image[0].size();

    // Apply the 1D Daubechies transform on each row
    for (int i = 0; i < rows; i++) {
        daubechies1D(image[i]);
    }
    // Apply the 1D Daubechies transform on each column
    for (int j = 0; j < cols; j++) {
        vector<float> column(rows);
        for (int i = 0; i < rows; i++) {
            column[i] = image[i][j];
        }
        daubechies1D(column);
        for (int i = 0; i < rows; i++) {
            image[i][j] = column[i];
        }
    }
}
// ***********************************************************************************************
// ********************************** Thresholding ***********************************************
void thresholdWaveletCoefficients(vector<vector<float>>& image, float threshold) {
    int rows = image.size();
    int cols = image[0].size();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Apply hard thresholding: Set coefficients below threshold to 0
            if (abs(image[i][j]) < threshold) {
                image[i][j] = 0.0f;
            }
        }
    }
}
// ********************************** Thresholding ***********************************************
// ***********************************************************************************************

// Function to apply wavelet compression
void waveletTransform(vector<u8>& buffer, u32 width, u32 height, double sigma) {
    vector<double> tempBuffer(width * height);

    // Convert the u8 buffer to double for the wavelet transform
    for (u32 i = 0; i < width * height; i++) {
        tempBuffer[i] = static_cast<double>(buffer[i]);
    }

    // Create a 2D image matrix
    vector<vector<float>> image(height, vector<float>(width));
    for (u32 row = 0; row < height; row++) {
        for (u32 col = 0; col < width; col++) {
            image[row][col] = static_cast<float>(tempBuffer[row * width + col]);
        }
    }

    // Apply the 2D Daubechies transform
    daubechies2D(image);


    // ****************************************************************************************
    // ************************************* Thresholding *************************************
    // Apply thresholding to the wavelet coefficients
    float threshold = 10.0f;  // You can adjust the threshold value based on desired compression
    thresholdWaveletCoefficients(image, threshold);
    // ************************************* Thresholding *************************************
    // ****************************************************************************************


    // Convert back to 1D buffer
    for (u32 row = 0; row < height; row++) {
        for (u32 col = 0; col < width; col++) {
            tempBuffer[row * width + col] = image[row][col];
        }
    }

    // Normalize to 0-255 range for output (8-bit)
    double minVal = *min_element(tempBuffer.begin(), tempBuffer.end());
    double maxVal = *max_element(tempBuffer.begin(), tempBuffer.end());

    for (u32 i = 0; i < width * height; i++) {
        double normalized = (tempBuffer[i] - minVal) / (maxVal - minVal) * 255.0;
        buffer[i] = static_cast<u8>(std::clamp(normalized, 0.0, 255.0));
    }
}

// Function to write the TIFF image from the buffer
void writeTiffImage(const char* filename, const vector<u8>& buffer, u32 width, u32 height, u16 samples_per_pixel) {
    TIFF* tiff = TIFFOpen(filename, "w");
    if (tiff == nullptr) {
        cerr << "Could not open TIFF file for writing!" << endl;
        return;
    }

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);  
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);  
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);

    for (u32 row = 0; row < height; row++) {
        TIFFWriteScanline(tiff, const_cast<u8*>(&buffer[row * width]), row);
    }

    TIFFClose(tiff);
}

void inverseDaubechies1D(vector<float>& data) {
    int N = data.size();
    vector<float> temp(N);

    for (int i = 0; i < N / 2; i++) {
        temp[2 * i] = ih0 * data[i] + ih1 * data[i + N / 2];
        temp[2 * i + 1] = ih2 * data[i] + ih3 * data[i + N / 2];
    }

    for (int i = 0; i < N; i++) {
        data[i] = temp[i];
    }
}

void inverseDaubechies2D(vector<vector<float>>& image) {
    int rows = image.size();
    int cols = image[0].size();

    // Apply the inverse 1D Daubechies transform on each row
    for (int i = 0; i < rows; i++) {
        inverseDaubechies1D(image[i]);
    }

    // Apply the inverse 1D Daubechies transform on each column
    for (int j = 0; j < cols; j++) {
        vector<float> column(rows);
        for (int i = 0; i < rows; i++) {
            column[i] = image[i][j];
        }
        inverseDaubechies1D(column);
        for (int i = 0; i < rows; i++) {
            image[i][j] = column[i];
        }
    }
}

void inverseWaveletTransform(vector<u8>& buffer, u32 width, u32 height) {
    // Convert the u8 buffer to double for inverse transform
    vector<double> tempBuffer(width * height);

    for (u32 i = 0; i < width * height; i++) {
        tempBuffer[i] = static_cast<double>(buffer[i]);
    }

    // Create a 2D image matrix for inverse transform
    vector<vector<float>> image(height, vector<float>(width));
    for (u32 row = 0; row < height; row++) {
        for (u32 col = 0; col < width; col++) {
            image[row][col] = static_cast<float>(tempBuffer[row * width + col]);
        }
    }

    // Apply the inverse 2D Daubechies transform
    inverseDaubechies2D(image);

    // Convert back to 1D buffer after the inverse transform
    for (u32 row = 0; row < height; row++) {
        for (u32 col = 0; col < width; col++) {
            tempBuffer[row * width + col] = image[row][col];
        }
    }

    // Normalize to 0-255 range for output (8-bit)
    double minVal = *min_element(tempBuffer.begin(), tempBuffer.end());
    double maxVal = *max_element(tempBuffer.begin(), tempBuffer.end());

    for (u32 i = 0; i < width * height; i++) {
        double normalized = (tempBuffer[i] - minVal) / (maxVal - minVal) * 255.0;
        buffer[i] = static_cast<u8>(std::clamp(normalized, 0.0, 255.0));
    }
}

int main() {

// ***********************
// Summary of the Process:
// ***********************
// 1) Read the TIFF Image: The image is read into a buffer (1D vector) and its dimensions are retrieved.
// 2) Apply Wavelet Transform (Data Transformation):
//    2a) Convert to Numeric Data: The 1D buffer is converted to a 2D array (image) of float values for the wavelet transform.
//    2b) Apply Wavelet Transform: The 2D Daubechies wavelet transform is applied in both row and column directions.
// 3) Write the Transformed TIFF File: The transformed image data is saved to a new TIFF file (using LZW compression).
// 4) Post-process and Normalize: The transformed data is normalized back to the 0-255 range and saved as an 8-bit TIFF image.
//    (This step also includes inverse wavelet transformation for reconstruction, which is optional but demonstrated.)


    const char* inputFile = "cameraman.tif";
    const char* outputFile = "output_compressed.tif";

    u32 width, height;
    u16 samples_per_pixel = 1;
    double sigma = 42; // Example sigma, change as needed
    vector<u8> buffer;

    // Step 1: Read the TIFF image
    readTiffImage(inputFile, buffer, width, height, samples_per_pixel);

    // Step 2: Apply wavelet compression
    waveletTransform(buffer, width, height, sigma);

    // Step 3: Write the compressed grayscale image to a new TIFF file
    writeTiffImage(outputFile, buffer, width, height, samples_per_pixel);

    // Step 4: Apply the inverse wavelet transform to reconstruct the image
    inverseWaveletTransform(buffer, width, height);

    // Optional: Save the reconstructed image to a new file
    const char* outputReconstructedFile = "output_reconstructed.tif";
    writeTiffImage(outputReconstructedFile, buffer, width, height, samples_per_pixel);

    cout << "Wavelet compression and reconstruction completed!" << endl;
    return 0;
}
