import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TfGroupAccordionComponent } from './tf-group-accordion.component';

describe('TfGroupAccordionComponent', () => {
  let component: TfGroupAccordionComponent;
  let fixture: ComponentFixture<TfGroupAccordionComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ TfGroupAccordionComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(TfGroupAccordionComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
