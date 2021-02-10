import Slider from '@material-ui/core/Slider';
import React, { FunctionComponent, useCallback } from 'react';

function valuetext(value: number) {
  return `Slice ${value}`;
}

type Props = {
  width: number
  numSlices: number
  currentSlice: number
  onCurrentSliceChanged: (x: number) => void
}

const SliceSlider: FunctionComponent<Props> = ({ width, numSlices, currentSlice, onCurrentSliceChanged }) => {

  const handleChange = useCallback(
    (event: React.ChangeEvent<{}>, value: number | number[]) => {
      onCurrentSliceChanged(value as any as number)
    },
    [onCurrentSliceChanged],
  )

  return (
    <div style={{ width }}>
      <Slider
        value={currentSlice}
        onChange={handleChange}
        getAriaValueText={valuetext}
        valueLabelDisplay="auto"
        step={1}
        marks
        min={0}
        max={numSlices - 1}
      />
    </div>
  );
}

export default SliceSlider