#include <iostream>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/video/video.hpp>

#include <vision/data/util.h>
#include <vision/srgb.h>

using namespace std;

int main_program(int argc, char** argv) {
	using namespace cv;

	Mat image = imread(vision::data::get_directory() + "/gamma_dalai_lama_gray.jpg", IMREAD_COLOR);
	imshow("Original", image);

	Mat smaller;
	resize(image, smaller, Size(), 0.5, 0.5, INTER_LANCZOS4);
	imshow("OpenCV Resize", smaller);

	Mat linear = vision::srgb_to_linear(image);
	Mat smaller_linear;
	resize(linear, smaller_linear, Size(), 0.5, 0.5, INTER_LANCZOS4);
	auto srgb_smaller_linear = vision::linear_to_srgb(smaller_linear);
	imshow("Petter Resize", srgb_smaller_linear);

	waitKey(0);
	return 0;
}

int main(int num_args, char* args[]) {
	try {
		auto res = main_program(num_args, args);
		cv::destroyAllWindows();
		return res;
	} catch (std::exception& err) {
		cerr << "Error: " << err.what() << endl;
		return 1;
	}
}
