/*
 * Copyright 2011, Blender Foundation.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor:
 *      Jeroen Bakker
 *      Monique Dewanchand
 */

#include "COM_BrightnessNode.h"
#include "COM_BrightnessOperation.h"
#include "COM_ExecutionSystem.h"

BrightnessNode::BrightnessNode(bNode *editorNode) : Node(editorNode)
{
	/* pass */
}

void BrightnessNode::convertToOperations(NodeConverter &converter, const CompositorContext &/*context*/) const
{
	bNode *bnode = this->getbNode();
	BrightnessOperation *operation = new BrightnessOperation();
	operation->setUsePremultiply((bnode->custom1 & 1) != 0);
	converter.addOperation(operation);

	converter.mapInputSocket(getInputSocket(0), operation->getInputSocket(0));
	converter.mapInputSocket(getInputSocket(1), operation->getInputSocket(1));
	converter.mapInputSocket(getInputSocket(2), operation->getInputSocket(2));
	converter.mapOutputSocket(getOutputSocket(0), operation->getOutputSocket(0));
}
