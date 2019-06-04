import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route, HashRouter } from "react-router-dom";
import Header from './Header.jsx';
import Navigation from './Navigation.jsx'
import Footer from './Footer.jsx';
import upload_contigs from './UploadContigs.jsx';
import Profile_view from './ProfileView.jsx';
import Upload_profile from './UploadProfile.jsx';
import Dendrogram_view from './DendrogramView.jsx';
import Tracking from './Tracking.jsx';
import Tutorial from './Tutorial.jsx';
import Tracking_result from './TrackingResult.jsx';
import About from './About.jsx';

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
                    <Route path="/cgMLST" exact component={upload_contigs} />
                    <Route path="/cgMLST/profiling" exact component={upload_contigs} />
                    <Route path="/cgMLST/profile_result" component={Profile_view} />
                    <Route path="/cgMLST/clustering" component={Upload_profile} />
                    <Route path="/cgMLST/clustering_result" component={Dendrogram_view} />
                    <Route path="/cgMLST/tracking" component={Tracking} />
                    <Route path="/cgMLST/tracking_result" component={Tracking_result} />
                    <Route path="/cgMLST/about" component={About} />
                    <Route path="/cgMLST/tutorial" component={Tutorial} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main />, document.getElementById('root'));
