import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from "react-router-dom";
import Header from './Header.jsx';
import Navigation from './Navigation.jsx'
import Footer from './Footer.jsx';
import upload_contigs from './UploadContigs.jsx';
import Profile_view from './ProfileView.jsx';
import Upload_profile from './UploadProfile.jsx';
import Dendrogram_view from './DendrogramView.jsx';
import Tracking from './Tracking.jsx';
import Example from './Example.jsx';
import Tutorial from './Tutorial.jsx';
import Tracking_result from './TrackingResult.jsx';
import Tracking_search from './TrackingSearch.jsx';
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
                    <Route path="/" exact component={upload_contigs} />
                    <Route path="/profile_view" component={Profile_view} />
                    <Route path="/upload_profile" component={Upload_profile} />
                    <Route path="/dendrogram_view" component={Dendrogram_view} />
                    <Route path="/tracking" component={Tracking} />
                    <Route path="/tracking_result" component={Tracking_result} />
                    <Route path="/tracking_search" component={Tracking_search} />
                    <Route path="/about" exact component={About} />
                    <Route path="/demo" component={Example} />
                    <Route path="/tutorial" component={Tutorial} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main /> ,document.getElementById('root'));
