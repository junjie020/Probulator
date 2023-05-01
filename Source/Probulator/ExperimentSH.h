#pragma once

#include <Probulator/Experiments.h>

namespace Probulator
{
template <size_t L> class ExperimentSH : public Experiment
{
public:
	void run(SharedData& data) override
	{
		SphericalHarmonicsT<vec3, L> shRadiance = {};

		const ivec2 imageSize = data.m_outputSize;
		data.m_directionImage.forPixels2D([&](const vec3& direction, ivec2 pixelPos)
		{
			float texelArea = latLongTexelArea(pixelPos, imageSize);
			vec3 radiance = (vec3)data.m_radianceImage.at(pixelPos);
			shAddWeighted(shRadiance, shEvaluate<L>(direction), radiance * texelArea);
		});

		if (m_targetLaplacian > 0.0f)
		{
			SphericalHarmonicsT<float, L> sh = shLuminance(shRadiance);
			m_lambda = shFindWindowingLambda(sh, m_targetLaplacian);
		}

		if (m_lambda != 0.0f)
		{
			shApplyWindowing<vec3, L>(shRadiance, m_lambda);
		}

		m_radianceImage = Image(data.m_outputSize);
		m_irradianceImage = Image(data.m_outputSize);

		// const float A[5] = {
		// 	pi,
		// 	pi * 2.0f / 3.0f,
		// 	pi * 1.0f / 4.0f,
		// 	0.0f,
		// 	-pi * 1.0f / 24.0f
		// };
		// for (int l = 0; l < L; ++l)
		// {
		// 	switch(l){
		// 	case 4:
		// 		for (int ii=0; ii<9; ++ii){
		// 			m_shCoeffs[16+ii] = vec4(shRadiance.data[15+ii] * A[4], 0.f);
		// 		}
		// 		break;
		// 	case 3:
		// 		for (int ii=0; ii<7; ++ii){
		// 			m_shCoeffs[9+ii] = vec4(shRadiance.data[9+ii] * A[3], 0.f);
		// 		}
		// 		break;
		// 	case 2:
		// 		for (int ii=0; ii<5; ++ii){
		// 			m_shCoeffs[4+ii] = vec4(shRadiance.data[4+ii] * A[2], 0.f);
		// 		}
		// 		break;
		// 	case 1:
		// 		for (int ii=0; ii<3; ++ii){
		// 			m_shCoeffs[1+ii] = vec4(shRadiance.data[1+ii] * A[1], 0.f);
		// 		}
		// 		break;
		// 	default:
		// 	case 0:
		// 		m_shCoeffs[0] = vec4(shRadiance.data[0] * A[0], 0.f); break;
		// 	}
		// }

		const float A[5] = {
			pi,
			pi * 2.0f / 3.0f,
			pi * 1.0f / 4.0f,
			0.0f,
			-pi * 1.0f / 24.0f
		};

		const int counts[5] = {1, 3, 5, 7, 9};
		int offset = 0;
		for (int l=0; l<L; ++l){
			int count = counts[l];
			for (int ii=offset; ii<count; ++ii){
				m_shCoeffs[ii] = vec4(shRadiance.data[ii] * A[l], 0.f);
			}
			offset += count;
		}

		// const int n = (L+1)*(L+1);
		// for (int ii=0; ii<n; ++ii){
		// 	m_shCoeffs[ii] = vec4(shRadiance.data[ii], 0.f);
		// }

		// auto irradiance_coeffs = [&shRadiance](vec3 n)
		// {
		// 	float c1 = 0.429043f;
		// 	float c2 = 0.511664f;
		// 	float c3 = 0.743125f;
		// 	float c4 = 0.886227f;
		// 	float c5 = 0.247708f;

		// 	vec3 L00;
		// 	vec3 L1_1, L10, L11;
		// 	vec3 L2_2, L2_1, L20, L21, L22;
		// 	//get_lighting_SH(L00, L1_1, L10, L11, L2_2, L2_1, L20, L21, L22);
		// 	L00  = shRadiance.data[0];

		// 	L1_1 = shRadiance.data[1];
		// 	L10  = shRadiance.data[2];
		// 	L11  = shRadiance.data[3];

		// 	L2_2 = shRadiance.data[4];
		// 	L2_1 = shRadiance.data[5];
		// 	L20  = shRadiance.data[6];
		// 	L21  = shRadiance.data[7];
		// 	L22  = shRadiance.data[8];

		// 	float x = n.x, y = n.y, z = n.z;
		// 	float x2 = x*x, y2 = y*y, z2 = z*z;
		// 	float xy = x*y, yz = y*z, xz = x*z;

		// 	return c1*L22*(x2-y2) + c3*L20*z2 + c4*L00 - c5*L20 
		// 			+ 2*c1*(L2_2*xy + L21*xz + L2_1*yz) 
		// 			+ 2*c2*(L11*x+L1_1*y+L10*z) ;
		// };

		// auto irradiance_coeffs = [this](vec3 n)
		// {
		// 	vec3 L00;
		// 	vec3 L1_1, L10, L11;
		// 	vec3 L2_2, L2_1, L20, L21, L22;
		// 	//get_lighting_SH(L00, L1_1, L10, L11, L2_2, L2_1, L20, L21, L22);
		// 	L00  = vec3(m_shCoeffs[0]);

		// 	L1_1 = vec3(m_shCoeffs[1]);
		// 	L10  = vec3(m_shCoeffs[2]);
		// 	L11  = vec3(m_shCoeffs[3]);

		// 	L2_2 = vec3(m_shCoeffs[4]);
		// 	L2_1 = vec3(m_shCoeffs[5]);
		// 	L20  = vec3(m_shCoeffs[6]);
		// 	L21  = vec3(m_shCoeffs[7]);
		// 	L22  = vec3(m_shCoeffs[8]);

		// 	const float x = -n.x;
		// 	const float y = -n.y;
		// 	const float z =  n.z;
		// 	float x2 = x*x, y2 = y*y, z2 = z*z;

		// 	vec3 result = L00 * 1.0f/(2.0f*sqrt(pi));;

		// 	result +=	L1_1 * -sqrt(3.0f/(4.0f*pi))*y + 
		// 				L10  *  sqrt(3.0f/(4.0f*pi))*z + 
		// 				L11  * -sqrt(3.0f/(4.0f*pi))*x;

		// 	result +=	L2_2 *  sqrt(15.0f/( 4.0f*pi))*y*x +
		// 				L2_1 * -sqrt(15.0f/( 4.0f*pi))*y*z +
		// 				L20  *  sqrt( 5.0f/(16.0f*pi))*(3.0f*z2-1.0f) +
		// 				L21  * -sqrt(15.0f/( 4.0f*pi))*x*z +
		// 				L22  *  sqrt(15.0f/(16.0f*pi))*(x2-y2);

		// 	return result;
		// };

		data.m_directionImage.forPixels2D([&](const vec3& direction, ivec2 pixelPos)
		{
			SphericalHarmonicsT<float, L> shDirection = shEvaluate<L>(direction);

			vec3 sampleSh = max(vec3(0.0f), shDot(shRadiance, shDirection));
			m_radianceImage.at(pixelPos) = vec4(sampleSh, 1.0f);

			vec3 sampleIrradianceSh = max(vec3(0.0f), shEvaluateDiffuse<vec3, L>(shRadiance, direction) / pi);

			//vec3 samplesh = max(vec3(0.f), irradiance_coeffs(direction) / pi);
			m_irradianceImage.at(pixelPos) = vec4(sampleIrradianceSh, 1.0f);
		});
	}

	void getProperties(std::vector<Property>& outProperties) override
	{
		Experiment::getProperties(outProperties);
		outProperties.push_back(Property("Lambda", &m_lambda));
		outProperties.push_back(Property("Target Laplacian", &m_targetLaplacian));
	}

	ExperimentSH<L>& setLambda(float v)
	{
		m_lambda = v;
		return *this;
	}

	ExperimentSH<L>& setTargetLaplacian(float v)
	{
		m_targetLaplacian = v;
		return *this;
	}

	float m_lambda = 0.0f;
	float m_targetLaplacian = -1.0f;
};

class ExperimentSHL1Geomerics : public Experiment
{
public:
	void run(SharedData& data) override
	{
		SphericalHarmonicsL1RGB shRadiance = {};

		const ivec2 imageSize = data.m_outputSize;
		data.m_directionImage.forPixels2D([&](const vec3& direction, ivec2 pixelPos)
		{
			float texelArea = latLongTexelArea(pixelPos, imageSize);
			vec3 radiance = (vec3)data.m_radianceImage.at(pixelPos);
			shAddWeighted(shRadiance, shEvaluateL1(direction), radiance * texelArea);
		});

		m_radianceImage = Image(data.m_outputSize);
		m_irradianceImage = Image(data.m_outputSize);

		data.m_directionImage.forPixels2D([&](const vec3& direction, ivec2 pixelPos)
		{
			SphericalHarmonicsL1 directionSh = shEvaluateL1(direction);

			vec3 sampleSh = max(vec3(0.0f), shDot(shRadiance, directionSh));
			m_radianceImage.at(pixelPos) = vec4(sampleSh, 1.0f);

			vec3 sampleIrradianceSh;
			for (u32 i = 0; i < 3; ++i)
			{
				SphericalHarmonicsL1 shRadianceChannel;
				shRadianceChannel[0] = shRadiance[0][i];
				shRadianceChannel[1] = shRadiance[1][i];
				shRadianceChannel[2] = shRadiance[2][i];
				shRadianceChannel[3] = shRadiance[3][i];
				sampleIrradianceSh[i] = shEvaluateDiffuseL1Geomerics(shRadianceChannel, direction) / pi;
			}
			m_irradianceImage.at(pixelPos) = vec4(sampleIrradianceSh, 1.0f);
		});
	}
};
}
