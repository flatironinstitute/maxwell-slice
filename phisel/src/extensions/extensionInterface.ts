import React from 'react';
import { CalculationPool } from "./common/hither";

export interface PhiselPlugin {
    name: string
    label: string
    priority?: number
    disabled?: boolean
    development?: boolean
    icon?: JSX.Element
}

export interface ViewPlugin extends PhiselPlugin {
    props?: {[key: string]: any}
    fullWidth?: boolean
    defaultExpanded?: boolean
    singleton?: boolean
}

interface ViewProps {
    calculationPool: CalculationPool
    width?: number
    height?: number
}

type SliceData = number[][]

type SlicedVolume = {
    dims: [number, number, number]
    res: [number, number, number]
    slices: SliceData[]
}

type Position = [number, number, number]

type SlicedVolumeSelection = {
    currentSlice: number
}

type SetCurrentSliceSlicedVolumeSelectionAction = {
    type: 'setCurrentSlice',
    currentSlice: number
}

export type SlicedVolumeSelectionAction = SetCurrentSliceSlicedVolumeSelectionAction

export type SlicedVolumeSelectionDispatch = (a: SlicedVolumeSelectionAction) => void

export interface SlicedVolumeViewProps extends ViewProps {
    volume: SlicedVolume
    selection: SlicedVolumeSelection
    selectionDispatch: SlicedVolumeSelectionDispatch
}

export interface SlicedVolumeViewPlugin extends ViewPlugin {
    component: React.ComponentType<SlicedVolumeViewProps>
    notebookCellHeight?: number
}

export const slicedVolumeSelectionReducer = (s: SlicedVolumeSelection, a: SlicedVolumeSelectionAction): SlicedVolumeSelection => {
    if (a.type === 'setCurrentSlice') {
        return {
            ...s,
            currentSlice: a.currentSlice
        }
    }
    else return s
}

export interface ExtensionContext {
    registerSlicedVolumeView: (V: SlicedVolumeViewPlugin) => void
}

export type PhiselPlugins = {
    slicedVolumeViews: SlicedVolumeViewPlugin[]
}

export const PluginsContext = React.createContext<PhiselPlugins>({slicedVolumeViews: []})