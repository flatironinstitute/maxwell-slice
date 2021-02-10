import { SlicedVolumeViewPlugin } from "./extensions/extensionInterface"

class ExtensionContextImpl {
    #slicedVolumeViewPlugins: SlicedVolumeViewPlugin[] = []
    registerSlicedVolumeView(p: SlicedVolumeViewPlugin) {
        this.#slicedVolumeViewPlugins.push(p)
    }
    plugins() {
        return {
            slicedVolumeViews: this.#slicedVolumeViewPlugins
        }
    }
}

const createExtensionContext = () => {
    return new ExtensionContextImpl()
}

export default createExtensionContext