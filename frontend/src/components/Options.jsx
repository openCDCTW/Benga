import React from 'react';
import ReactDOM from 'react-dom';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';

const styles = theme => ({
  root: {
    display: 'flex',
    flexWrap: 'wrap',
  },
  formControl: {
    margin: theme.spacing.unit * 5,
    minWidth: 240,
  },
  selectEmpty: {
    marginTop: theme.spacing.unit * 2,
  },
});

class Options extends React.Component {

    constructor(props) {
        super(props);

        this.state = {
            database: '',
            labelWidth: 0,
        };

        this.handleChange = this.handleChange.bind(this);
    }
    
    handleChange(event){
        this.setState(state => ({ [event.target.name]: event.target.value}));
        window.databaseName = event.target.value;
    };

    render(){

        const { classes } = this.props;
        
        return (
          <form className={classes.root} autoComplete="off">
            <FormControl required className={classes.formControl} disabled={this.props.switch}>
              <InputLabel htmlFor="database-required">Database</InputLabel>
              <Select
                value={this.state.database}
                onChange={this.handleChange}
                name="database"
                inputProps={{
                  id: 'database-required',
                }}
                className={classes.selectEmpty}
                >
                <MenuItem value="">
                  <em>Choose Database</em>
                </MenuItem>
                <MenuItem value={'Salmonella_enterica'}>Salmonella_enterica</MenuItem>
                <MenuItem value={'Vibrio_cholerae'}>Vibrio cholerae</MenuItem>
                <MenuItem value={'Listeria_monocytogenes'}>Listeria monocytogenes</MenuItem>
                <MenuItem value={'Campylobacter_jejuni'}>Campylobacter jejuni</MenuItem>
              </Select>
              <FormHelperText>Required</FormHelperText>
            </FormControl>
          </form>
        );
    }
}

Options.propTypes = {
  classes: PropTypes.object.isRequired,
};

export default withStyles(styles)(Options);