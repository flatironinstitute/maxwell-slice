import React from 'react';
import ReactDOM from 'react-dom';
import App from './App';
import createExtensionContext from './createExtensionContext';
import { PluginsContext } from './extensions/extensionInterface';
import './index.css';
import registerExtensions from './registerExtensions';
import reportWebVitals from './reportWebVitals';

const extensionContext = createExtensionContext()
registerExtensions(extensionContext)

ReactDOM.render(
  <React.StrictMode>
    <PluginsContext.Provider value={extensionContext.plugins()}>
      <App />
    </PluginsContext.Provider>
  </React.StrictMode>,
  document.getElementById('root')
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();
