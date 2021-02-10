// PHISEL-EXTENSION: slicedvolumeview1
// PHISEL-EXTENSION-TAGS: jupyter

import { BubbleChart } from '@material-ui/icons';
import React from 'react';
import { ExtensionContext } from '../extensionInterface';
import SlicedVolumeView1 from './SlicedVolumeView1/SlicedVolumeView1';

export function activate(context: ExtensionContext) {
    context.registerSlicedVolumeView({
        name: 'SliceView',
        label: 'Slice view',
        priority: 50,
        component: SlicedVolumeView1,
        icon: <BubbleChart />
    })
}