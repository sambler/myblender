/*******************************************************************************
* Copyright 2009-2015 Juan Francisco Crespo Galán
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

#include "Exception.h"

#include <cassert>

using namespace aud;

#define AUD_CAPI_IMPLEMENTATION
#include "AUD_ThreadPool.h"

AUD_API AUD_ThreadPool* AUD_ThreadPool_create(int nThreads)
{
	try
	{
		return new AUD_ThreadPool(new ThreadPool(nThreads));
	}
	catch(Exception&)
	{
		return nullptr;
	}
}

AUD_API void AUD_ThreadPool_free(AUD_ThreadPool* pool)
{
	assert(pool);
	delete pool;
}