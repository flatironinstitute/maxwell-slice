import React, { FunctionComponent } from 'react';
import CanvasWidget from '../../common/CanvasWidget';
import { useLayer, useLayers } from '../../common/CanvasWidget/CanvasWidgetLayer';
import { createMainLayer } from './mainLayer';

type Props = {
    width: number
    height: number
    slice: number[][]
    valueRange: {min: number, max: number}
}

const SliceView: FunctionComponent<Props> = ({width, height, slice, valueRange}) => {
    
    const mainLayerProps = {
        width: 300,
        height: 300,
        slice: slice,
        valueRange: valueRange
    }
    const mainLayer = useLayer(createMainLayer, mainLayerProps)
    const layers = useLayers([mainLayer])

    return (
        <CanvasWidget
            layers={layers}
            width={300}
            height={300}
        />
    )
}

export default SliceView