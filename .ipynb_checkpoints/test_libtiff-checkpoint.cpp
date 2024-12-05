#include <tiffio.h>
#include <iostream>
#include <vector>

void readTiffImage(const char* filename) {
    TIFF* tiff = TIFFOpen(filename, "r");
    if (tiff == nullptr) {
        std::cerr << "Could not open TIFF file!" << std::endl;
        return;
    }

    // Get image width and height
    uint32_t width, height;
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);

    // Get the number of channels (assuming 8-bit per channel RGB)
    uint16_t samples_per_pixel;
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

    // Allocate memory for the image data
    std::vector<uint8_t> buffer(width * height * samples_per_pixel);

    // Read image into buffer
    for (uint32_t row = 0; row < height; row++) {
        TIFFReadScanline(tiff, &buffer[row * width * samples_per_pixel], row);
    }

    TIFFClose(tiff);

    std::cout << "Image read: " << width << "x" << height << " with " << samples_per_pixel << " channels." << std::endl;
    // Process the image data in the buffer (next steps).
}

int main() {
    const char* filename = "board.tif";
    readTiffImage(filename);
    return 0;
}