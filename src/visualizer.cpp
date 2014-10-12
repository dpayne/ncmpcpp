/***************************************************************************
 *   Copyright (C) 2008-2014 by Andrzej Rybczak                            *
 *   electricityispower@gmail.com                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.              *
 ***************************************************************************/

#include "visualizer.h"

#ifdef ENABLE_VISUALIZER

#include <boost/date_time/posix_time/posix_time.hpp>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <fcntl.h>

#include "global.h"
#include "settings.h"
#include "status.h"
#include "statusbar.h"
#include "title.h"
#include "screen_switcher.h"
#include "status.h"

using Global::MainStartY;
using Global::MainHeight;

Visualizer *myVisualizer;

namespace {

const int fps = 25;

const Color colorMap [] = { clWhite, clCyan, clBlue, clGreen, clYellow, clMagenta, clRed };
const char asciiGreyScaleMap [] = { '@', '%', '#', '*', '+', '=', '-' };

}

Visualizer::Visualizer()
: Screen(NC::Window(0, MainStartY, COLS, MainHeight, "", Config.visualizer_color, NC::Border::None))
{
	ResetFD();
	m_samples = 44100/fps;
	if (Config.visualizer_in_stereo)
		m_samples *= 2;
#	ifdef HAVE_FFTW3_H
	m_fftw_results = m_samples/2+1;
	m_freq_magnitudes = new double[m_fftw_results];
	m_fftw_input = static_cast<double *>(fftw_malloc(sizeof(double)*m_samples));
	m_fftw_output = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*m_fftw_results));
	m_fftw_plan = fftw_plan_dft_r2c_1d(m_samples, m_fftw_input, m_fftw_output, FFTW_ESTIMATE);
#	endif // HAVE_FFTW3_H
}

void Visualizer::switchTo()
{
	SwitchTo::execute(this);
	w.clear();
	SetFD();
	m_timer = boost::posix_time::from_time_t(0);
	drawHeader();
}

void Visualizer::resize()
{
	size_t x_offset, width;
	getWindowResizeParams(x_offset, width);
	w.resize(width, MainHeight);
	w.moveTo(x_offset, MainStartY);
	hasToBeResized = 0;
}

std::wstring Visualizer::title()
{
	return L"Music visualizer";
}

void Visualizer::update()
{
	if (m_fifo < 0)
		return;

	// PCM in format 44100:16:1 (for mono visualization) and 44100:16:2 (for stereo visualization) is supported
	int16_t buf[m_samples];
	ssize_t data = read(m_fifo, buf, sizeof(buf));
	if (data < 0) // no data available in fifo
		return;

	if (m_output_id != -1 && Global::Timer - m_timer > Config.visualizer_sync_interval)
	{
		Mpd.DisableOutput(m_output_id);
		usleep(50000);
		Mpd.EnableOutput(m_output_id);
		m_timer = Global::Timer;
	}

	void (Visualizer::*draw)(int16_t *, ssize_t, size_t, size_t);
#	ifndef HAVE_FFTW3_H
	draw = &Visualizer::DrawSoundWave;
#	else
	if ( Config.visualizer_function == 1 )
	{
		draw = &Visualizer::DrawSoundWaveFillAscii;
		color = false;
	}
	else if ( Config.visualizer_function == 2 )
	{
		draw = &Visualizer::DrawSoundWaveFill;
		color = true;
	}
	else if ( Config.visualizer_function == 3 )
	{
		draw = &Visualizer::DrawSoundWaveFillAscii;
		color = true;
	}
	else if ( Config.visualizer_function == 4 )
	{
		draw = &Visualizer::DrawFrequencySpectrum;
		color = false;
	}
	else if ( Config.visualizer_function == 5 )
	{
		draw = &Visualizer::DrawFrequencySpectrum;
		color = true;
	}
	else
	{
		draw = &Visualizer::DrawSoundWave;
		color = false;
	}
#	endif // HAVE_FFTW3_H

	const ssize_t samples_read = data/sizeof(int16_t);
	std::for_each(buf, buf+samples_read, [](int16_t &sample) {
		int32_t tmp = sample * Config.visualizer_sample_multiplier;
		if (tmp < std::numeric_limits<int16_t>::min())
			sample = std::numeric_limits<int16_t>::min();
		else if (tmp > std::numeric_limits<int16_t>::max())
			sample = std::numeric_limits<int16_t>::max();
		else
			sample = tmp;
	});

	w.clear();
	if (Config.visualizer_in_stereo)
	{
		int16_t buf_left[samples_read/2], buf_right[samples_read/2];
		for (ssize_t i = 0, j = 0; i < samples_read; i += 2, ++j)
		{
			buf_left[j] = buf[i];
			buf_right[j] = buf[i+1];
		}
		size_t half_height = MainHeight/2;
		(this->*draw)(buf_left, samples_read/2, 0, half_height);
		(this->*draw)(buf_right, samples_read/2, half_height+(draw == &Visualizer::DrawSoundWave ? 1 : 0), half_height+(draw != &Visualizer::DrawSoundWave ? 1 : 0));
	}
	else
		(this->*draw)(buf, samples_read, 0, MainHeight);
	w.refresh();
}

int Visualizer::windowTimeout()
{
	if (m_fifo >= 0 && Status::State::player() == MPD::psPlay)
		return 1000/fps;
	else
		return Screen<WindowType>::windowTimeout();
}

void Visualizer::spacePressed()
{
#	ifdef HAVE_FFTW3_H
	Config.visualizer_use_wave = !Config.visualizer_use_wave;
	Statusbar::printf("Visualization type: %1%",
		Config.visualizer_use_wave ? "sound wave" : "frequency spectrum"
	);
#	endif // HAVE_FFTW3_H
}

Color Visualizer::toColor( int number, int max )
{
	const int normalizedNumber = ( ( number * 7 ) / max ) % 7;
	return colorMap[ normalizedNumber ];
}

char Visualizer::toAsciiGrey( int number, int max )
{
	const int normalizedNumber = ( ( number * 7 ) / max ) % 7;
	return asciiGreyScaleMap[ normalizedNumber ];
}

void Visualizer::DrawSoundWave(int16_t *buf, ssize_t samples, size_t y_offset, size_t height, bool color)
{
	const int samples_per_col = samples/w->GetWidth();
	const int half_height = height/2;
	double prev_point_pos = 0;
	const size_t win_width = w->GetWidth();
	for (size_t i = 0; i < win_width; ++i)
	{
		double point_pos = 0;
		for (int j = 0; j < samples_per_col; ++j)
			point_pos += buf[i*samples_per_col+j];
		point_pos /= samples_per_col;
		point_pos /= std::numeric_limits<int16_t>::max();
		point_pos *= half_height;
		*w << NC::XY(i, y_offset+half_height+point_pos) << Config.visualizer_chars[0];
		if (i && abs(prev_point_pos-point_pos) > 2)
		{
			// if gap is too big. intermediate values are needed
			// since without them all we see are blinking points
			const int breakpoint = std::max(prev_point_pos, point_pos);
			const int half = (prev_point_pos+point_pos)/2;
			*w << clDefault;
			for (int k = std::min(prev_point_pos, point_pos)+1; k < breakpoint; k += 2)
				*w << NC::XY(i-(k < half), y_offset+half_height+k) << Config.visualizer_chars[0];
		}
		prev_point_pos = point_pos;
	}
}

void Visualizer::DrawSoundWaveFill(int16_t *buf, ssize_t samples, size_t y_offset, size_t height, bool color)
{
	const int samples_per_col = samples/w->GetWidth();
	const int half_height = height/2;
	double prev_point_pos = 0;
	const size_t win_width = w->GetWidth();
	const bool left = y_offset > 0;
	int x = 0;
	for (size_t i = 0; i < win_width; ++i)
	{
		double point_pos = 0;
		for (int j = 0; j < samples_per_col; ++j)
			point_pos += buf[i*samples_per_col+j];
		point_pos /= samples_per_col;
		point_pos /= std::numeric_limits<int16_t>::max();
		point_pos *= half_height;
		for (int k = 0; k < point_pos * 2; k += 2)
		{
			x = height;
			if ( left )
			{
				x += k;
			}
			else
			{
				x -= k;
			}

			if ( x > 0 && x < w->GetHeight() && (i-(k < half_height + point_pos)) > 0 && (i-(k < half_height + point_pos)) < w->GetWidth() )
			{
				if ( color )
				{
					*w << toColor( k, height );
				}
				else
				{
					*w << clDefault;
				}
				*w << NC::XY(i-(k < half_height + point_pos), x) << Config.visualizer_chars[0];
			}
		}
		prev_point_pos = point_pos;
	}
}

void Visualizer::DrawSoundWaveFillAscii(int16_t *buf, ssize_t samples, size_t y_offset, size_t height, bool color)
{
	const int samples_per_col = samples/w->GetWidth();
	const int half_height = height/2;
	double prev_point_pos = 0;
	const size_t win_width = w->GetWidth();
	const bool left = y_offset > 0;
	int x = 0;
	for (size_t i = 0; i < win_width; ++i)
	{
		double point_pos = 0;
		for (int j = 0; j < samples_per_col; ++j)
			point_pos += buf[i*samples_per_col+j];
		point_pos /= samples_per_col;
		point_pos /= std::numeric_limits<int16_t>::max();
		point_pos *= half_height;
		for (int k = 0; k < point_pos * 2; k += 2)
		{
			x = height;
			if ( left )
			{
				x += k;
			}
			else
			{
				x -= k;
			}

			if ( x > 0 && x < w->GetHeight() && (i-(k < half_height + point_pos)) > 0 && (i-(k < half_height + point_pos)) < w->GetWidth() )
			{
				if ( color )
				{
					*w << toColor( k, height );
				}
				else
				{
					*w << clDefault;
				}
				*w << NC::XY(i-(k < half_height + point_pos), x) << toAsciiGrey( k, height );
			}
		}
		prev_point_pos = point_pos;
	}
}

#ifdef HAVE_FFTW3_H
void Visualizer::DrawFrequencySpectrum(int16_t *buf, ssize_t samples, size_t y_offset, size_t height, bool color)
{
	for (unsigned i = 0, j = 0; i < itsSamples; ++i)
	{
		if (j < samples)
			itsInput[i] = buf[j++];
		else
			itsInput[i] = 0;
	}

	fftw_execute(itsPlan);

	// count magnitude of each frequency and scale it to fit the screen
	for (unsigned i = 0; i < itsFFTResults; ++i)
		itsFreqsMagnitude[i] = sqrt(itsOutput[i][0]*itsOutput[i][0] + itsOutput[i][1]*itsOutput[i][1])/1e5*height/5;

	const size_t win_width = w->GetWidth();
	const int freqs_per_col = (itsFFTResults/win_width) /* cut bandwidth a little to achieve better look */ * 7/10;
	for (size_t i = 0; i < win_width; i += 1)
	{
		size_t bar_height = 0;
		for (int j = 0; j < freqs_per_col; ++j)
			bar_height += itsFreqsMagnitude[i*freqs_per_col+j];

		bar_height = std::min(bar_height/freqs_per_col, height);
		const size_t start_y = y_offset > 0 ? y_offset : height-bar_height;
		const size_t stop_y = std::min(bar_height+start_y, w->GetHeight());
		const Color colorHeight = toColor( i, win_width );

		for (size_t j = start_y; j < stop_y; j += 1)
		{
			if ( color )
			{
				*w << colorHeight;
			}
			else
			{
				*w << clDefault;
			}
			*w << NC::XY(i, j) << Config.visualizer_chars[1];
		}
	}
}
#endif // HAVE_FFTW3_H

void Visualizer::SetFD()
{
	if (m_fifo < 0 && (m_fifo = open(Config.visualizer_fifo_path.c_str(), O_RDONLY | O_NONBLOCK)) < 0)
		Statusbar::printf("Couldn't open \"%1%\" for reading PCM data: %2%",
			Config.visualizer_fifo_path, strerror(errno)
		);
}

void Visualizer::ResetFD()
{
	m_fifo = -1;
}

void Visualizer::FindOutputID()
{
	m_output_id = -1;
	if (!Config.visualizer_output_name.empty())
	{
		size_t idx = 0;
		Mpd.GetOutputs([this, &idx](MPD::Output output) {
			if (output.name() == Config.visualizer_output_name)
				m_output_id = idx;
			++idx;
		});
		if (m_output_id == -1)
			Statusbar::printf("There is no output named \"%s\"", Config.visualizer_output_name);
	}
}

#endif // ENABLE_VISUALIZER

/* vim: set tabstop=4 softtabstop=4 shiftwidth=4 noexpandtab : */
