import React from 'react';
import { withStyles } from '@material-ui/core/styles';
import AppBar from '@material-ui/core/AppBar';
import CssBaseline from '@material-ui/core/CssBaseline';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';
import Button from '@material-ui/core/Button';
import PageContent from './PageContent.jsx'

const styles = theme => ({
  root: {
    display: 'flex',
  },
  title:{
    color: 'white',
    flexGrow: 1,
  },
  loginbtn:{
    color: 'white',
    display: 'none',
  },
  appBar: {
    zIndex: 1201,
  },
});

class Main extends React.Component {

    constructor(props) {
        super(props);
        history.pushState(null, null, location.href);
        window.onpopstate = function(){
            history.go(1);
        };
    }

     render() {
        const { classes } = this.props;

        return(
            <div className={classes.root}>
                <CssBaseline />
                <AppBar position="fixed" className={classes.appBar}>
                    <Toolbar>
                        <Typography className={classes.title} variant="h6">
                            cgMLST@Taiwan
                        </Typography>
                        <Button className={classes.loginbtn}>Login</Button>
                    </Toolbar>
                </AppBar>
                <PageContent />
            </div>
        );
    }
}

export default withStyles(styles)(Main);