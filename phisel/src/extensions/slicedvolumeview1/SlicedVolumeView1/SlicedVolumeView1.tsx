import React, { FunctionComponent, useCallback, useMemo } from 'react';
import { SlicedVolumeViewProps } from '../../extensionInterface';
import SliceSlider from './SliceSlider';
import SliceView from './SliceView';

const SlicedVolumeView1: FunctionComponent<SlicedVolumeViewProps> = ({volume, selection, selectionDispatch}) => {
    const slices = volume.slices
    const valueRange = useMemo(() => {
        const range = {min: 0, max: 0}
        let first = true
        for (let slice of slices) {
            for (let i = 0; i < slice.length; i++) {
                for (let j = 0; j < slice.length; j++) {
                    const v = slice[i][j]
                    range.min = first ? v : Math.min(range.min, v)
                    range.max = first ? v : Math.max(range.max, v)
                    first = false
                }
            }
        }
        return range
    }, [slices])

    const handleCurrentSliceChanged = useCallback((currentSlice: number) => {
            selectionDispatch({type: 'setCurrentSlice', currentSlice})
        },
        [selectionDispatch]
    )

    const width = 300

    return (
        <div style={{margin: 30}}>
            <SliceView
                width={width}
                height={300}
                slice={slices[selection.currentSlice]}
                valueRange={valueRange}
            />
            <SliceSlider
                width={width}
                numSlices={slices.length}
                currentSlice={selection.currentSlice}
                onCurrentSliceChanged={handleCurrentSliceChanged}
            />
        </div>
    )
}

export default SlicedVolumeView1