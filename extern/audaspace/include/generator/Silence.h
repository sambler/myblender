/*******************************************************************************
 * Copyright 2009-2016 Jörg Müller
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

#pragma once

/**
 * @file Silence.h
 * @ingroup generator
 * The Silence class.
 */

#include "ISound.h"
#include "respec/Specification.h"

AUD_NAMESPACE_BEGIN

/**
 * This sound creates a reader that plays silence.
 */
class AUD_API Silence : public ISound
{
private:
	/**
	 * The target sample rate for output.
	 */
	const SampleRate m_sampleRate;

	// delete copy constructor and operator=
	Silence(const Silence&) = delete;
	Silence& operator=(const Silence&) = delete;

public:
	/**
	 * Creates a new silence sound.
	 * \param sampleRate The target sample rate for playback.
	 */
	Silence(SampleRate sampleRate = RATE_48000);

	virtual std::shared_ptr<IReader> createReader();
};

AUD_NAMESPACE_END
