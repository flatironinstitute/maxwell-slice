import React, { FunctionComponent, useContext, useMemo, useReducer } from 'react';
import { createCalculationPool } from './extensions/common/hither';
import { PluginsContext, slicedVolumeSelectionReducer, SlicedVolumeViewProps } from './extensions/extensionInterface';

const calculationPool = createCalculationPool({maxSimultaneous: 6})

type Slice = number[][]

const TestPage: FunctionComponent<{}> = () => {
  const plugins = useContext(PluginsContext)

  const [selection, selectionDispatch] = useReducer(slicedVolumeSelectionReducer, {currentSlice: 0})

  const n1 = 200
  const n2 = 200
  const numSlices = 12
  const slices = useMemo(() => {
    const ret: Slice[] = []
    for (let i = 0; i < numSlices; i++) {
        ret.push(createRandomSlice(n1, n2))
    }
    return ret
  }, [n1, n2, numSlices])

  const props: SlicedVolumeViewProps = {
    volume: {
      dims: [3, 3, 3],
      res: [1, 1, 1],
      slices
    },
    selection,
    selectionDispatch,
    calculationPool
  }

  return (
    <div>
      {
        plugins.slicedVolumeViews.map(p => (
          <div><p.component {...props} /></div>
        ))
      }
    </div>
  );
}

const createRandomSlice = (n1: number, n2: number) => {
    const x: Slice = []
    for (let i1 = 0; i1 < n1; i1++) {
        const y: number[] = []
        for (let i2 = 0; i2 < n2; i2 ++) {
            y.push(Math.random())
        }
        x.push(y)
    }
    return x
}

export default TestPage;
