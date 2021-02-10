
// !begin-code-generation!
import { ExtensionContext } from './extensions/extensionInterface'
import { activate as activateslicedvolumeview1 } from './extensions/slicedvolumeview1/slicedvolumeview1'

// !end-code-generation!

/*
Extensions are automatically detected and added to this file via code generation (see task configured in vscode)
They must be .tsx files with the following appearing at the top of the file
// PHISEL-EXTENSION: <name>
And they must include an activate() function
Use the following to also include the extension in the jupyterlab extension:
// PHISEL-EXTENSION-TAGS: jupyter
*/


const registerExtensions = (context: ExtensionContext) => {
    // !begin-code-generation!
    activateslicedvolumeview1(context)
    // !end-code-generation!
}

export default registerExtensions

// !note! This file involves code generation.
// !note! The following template file was used: ./src/registerExtensions.ts.j2
// !note! You may edit the generated file outside of the code-generation blocks.
// !note! Changes to the generated file will be updated in the template file.
// !note! If vscode automatically moves the code-generation block delimiters upon save, just manually move them back and re-save
