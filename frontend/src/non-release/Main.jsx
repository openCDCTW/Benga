import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route, HashRouter } from "react-router-dom";
import Header from './Header.jsx';
import Navigation from '../components/Navigation.jsx'
import Footer from '../components/Footer.jsx';
import upload_contigs from '../components/UploadContigs.jsx';

class Main extends React.Component {

    constructor(props) {

        super(props);
        
        history.pushState(null, null, location.href);
        window.onpopstate = function(){
            history.go(1);
        };
    }

    render() {
        return(
            <Router>
                <div>
                    <Header />
                    <Navigation />
                    <Route path="/cgMLST/non-release/" exact component={upload_contigs} />
                    <Route path="/cgMLST/non-release/profiling" exact component={upload_contigs} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main />, document.getElementById('root'));
