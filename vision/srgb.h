
#include <cmath>

#include <opencv2/core/core.hpp>

namespace vision {
// Convertes an sRGB value in the range [0, 255] to a linear value in
// the range [0, 1].
inline float srgb_to_linear(float srgb) {
	auto linear = srgb / 255.0f;
	if (linear <= 0.4045f) {
		linear = linear / 12.92f;
	} else {
		linear = std::powf((linear + 0.055f) / 1.055f, 2.4f);
	}
	return linear;
}

// Converts a linear value in the range [0, 1] to an sRGB value in
// the range [0, 255].
inline float linear_to_srgb(float linear) {
	float srgb;
	if (linear <= 0.0031308f) {
		srgb = linear * 12.92f;
	} else {
		srgb = 1.055f * std::powf(linear, 1.0f / 2.4f) - 0.055f;
	}
	return srgb * 255.f;
}

// Converts an OpenCV Mat with sRGB values in the range [0, 255] to
// a float (CV_32FC3) OpenCV Mat in the range [0, 1].
static cv::Mat srgb_to_linear(const cv::Mat& mat) {
	cv::Mat linear;
	mat.convertTo(linear, CV_32FC3);
	for (auto itr = linear.begin<cv::Vec<float, 3>>(); itr != linear.end<cv::Vec<float, 3>>();
	     ++itr) {
		for (int channel = 0; channel < 3; ++channel) {
			(*itr)[channel] = srgb_to_linear((*itr)[channel]);
		}
	}
	return linear;
}

// Converts an OpenCV Mat with linear values in the range [0, 1] to
// a 8-bit (CV_8UC3) OpenCV Mat in the range [0, 255].
static cv::Mat linear_to_srgb(const cv::Mat& linear) {
	cv::Mat srgb(linear.rows, linear.cols, CV_8UC3);
	for (int r = 0; r < linear.rows; ++r) {
		for (int c = 0; c < linear.cols; ++c) {
			auto& dest = srgb.at<cv::Vec<unsigned char, 3>>(r, c);
			auto& source = linear.at<cv::Vec<float, 3>>(r, c);
			for (int channel = 0; channel < 3; ++channel) {
				dest[channel] =
				    cv::saturate_cast<unsigned char>(linear_to_srgb(source[channel]) + 0.5f);
			}
		}
	}
	return srgb;
}
}  // namespace vision
